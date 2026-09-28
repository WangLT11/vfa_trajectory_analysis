options(stringsAsFactors = FALSE, scipen = 999)

# Run second. This script rebuilds the 24-month df=2, K=3 trajectory model
# from raw data and saves the shared object used by scripts 03 and 04.

input_file <- Sys.getenv("VFA_DATA_FILE", "data/vfa_body_composition.csv")
input_encoding <- Sys.getenv("VFA_INPUT_ENCODING", "UTF-8")
output_dir <- Sys.getenv("VFA_OUTPUT_DIR", "results")
local_lib <- Sys.getenv("VFA_R_LIBRARY", "")
if (nzchar(local_lib) && dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

packages <- c("data.table", "dplyr", "tidyr", "ggplot2", "survival", "survminer", "gridExtra", "lcmm", "MASS")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(survival)
  library(survminer); library(gridExtra); library(lcmm); library(MASS)
})
if (.Platform$OS.type == "windows") try(windowsFonts(Arial = windowsFont("Arial")), silent = TRUE)

out <- file.path(output_dir, "02_trajectory_cox")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20260814)

colors <- c("High-stable" = "#2A9D8F", "Progressive-loss" = "#C84B31", "Rebound-increasing" = "#3A6EA5")
safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
safe_id <- function(x) {
  x <- trimws(iconv(as.character(x), from = "UTF-8", to = "UTF-8", sub = ""))
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_; x
}
parse_date <- function(x) {
  z <- trimws(iconv(as.character(x), from = "UTF-8", to = "UTF-8", sub = ""))
  z[z %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  ans <- rep(as.Date(NA), length(z)); n <- suppressWarnings(as.numeric(z))
  excel <- which(!is.na(n) & n > 20000 & n < 60000)
  if (length(excel)) ans[excel] <- as.Date(n[excel], origin = "1899-12-30")
  remaining <- which(is.na(ans) & !is.na(z))
  for (fmt in c("%Y/%m/%d", "%Y-%m-%d", "%Y/%m/%d %H:%M:%S", "%Y-%m-%d %H:%M:%S", "%Y%m%d")) {
    if (!length(remaining)) break
    parsed <- as.Date(z[remaining], format = fmt); ok <- !is.na(parsed)
    ans[remaining[ok]] <- parsed[ok]; remaining <- remaining[!ok]
  }
  ans
}
first_col <- function(nms, candidates, required = TRUE) {
  z <- candidates[candidates %in% nms]
  if (!length(z)) {
    if (required) stop("Missing column: ", paste(candidates, collapse = " | "))
    return(NA_character_)
  }
  z[1]
}
audit_extreme <- function(x) {
  ans <- x
  for (i in which(is.finite(ans) & ans > 400)) {
    med <- median(ans[-i][is.finite(ans[-i])], na.rm = TRUE)
    if (is.finite(med) && med > 0 && ans[i] / med >= 10) ans[i] <- NA_real_
  }
  ans
}
collapse_stage <- function(x) case_when(x == 1 ~ "I", x %in% 2:4 ~ "II", x %in% 5:7 ~ "III", TRUE ~ NA_character_)
theme_pub <- function(size = 10.5) {
  theme_classic(base_family = "Arial", base_size = size) +
    theme(panel.grid = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black"), legend.position = "bottom",
          plot.background = element_rect(fill = "white", colour = NA))
}
pool_scalar <- function(q, u) {
  m <- length(q); qbar <- mean(q); ubar <- mean(u); b <- var(q)
  total <- ubar + (1 + 1 / m) * b
  df <- if (!is.finite(b) || b < 1e-14) Inf else (m - 1) * (1 + ubar / ((1 + 1 / m) * b))^2
  se <- sqrt(total); crit <- if (is.finite(df)) qt(.975, df) else qnorm(.975)
  p <- if (is.finite(df)) 2 * pt(-abs(qbar / se), df) else 2 * pnorm(-abs(qbar / se))
  c(estimate = qbar, lower = qbar - crit * se, upper = qbar + crit * se, p = p)
}

raw <- data.table::fread(input_file, encoding = input_encoding, data.table = FALSE, check.names = FALSE)
nms <- names(raw)
id_col <- first_col(nms, c(" ID", "ID")); age_col <- first_col(nms, c("Age (years)", "Age"))
surgery_col <- first_col(nms, c("Surgery_Date", "Surgery Date")); os_col <- first_col(nms, "Overall_Survival_Months")
death_col <- first_col(nms, "Death_Event"); sex_col <- first_col(nms, c("Gender (1=Male, 2=Female)", "Gender"))
bmi_col <- first_col(nms, "BMI")
stage_col <- first_col(nms, c("pStage(1=I,2=IIA,3=IIB,4=IIC,5=IIIA,6=IIIB,7=IIIC)", "pStage"))
chemo_col <- first_col(nms, c("Adjuvant_Chemotherapy?(0=No, 1=Yes)", "Adjuvant_Chemotherapy"))

clinical <- data.frame(
  ID = safe_id(raw[[id_col]]), age = safe_num(raw[[age_col]]), surgery_date = parse_date(raw[[surgery_col]]),
  os_months = safe_num(raw[[os_col]]), death_raw = safe_num(raw[[death_col]]),
  sex_code = safe_num(raw[[sex_col]]), bmi = safe_num(raw[[bmi_col]]),
  stage_code = safe_num(raw[[stage_col]]), chemo_num = safe_num(raw[[chemo_col]])
) %>%
  mutate(
    death_event = case_when(death_raw == 0 ~ 0L, death_raw == 1 ~ 1L, TRUE ~ NA_integer_),
    sex = factor(case_when(sex_code == 1 ~ "Male", sex_code == 2 ~ "Female", TRUE ~ NA_character_),
                 levels = c("Male", "Female")),
    stage = factor(collapse_stage(stage_code), levels = c("I", "II", "III")),
    chemo_num = ifelse(chemo_num %in% c(0, 1), chemo_num, NA_real_),
    chemo = factor(chemo_num, levels = c(0, 1), labels = c("No", "Yes"))
  ) %>%
  filter(!is.na(ID), !is.na(age), age < 85, !is.na(surgery_date), !is.na(os_months), os_months > 0)

visits <- c("Pre", "1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th")
availability <- bind_rows(lapply(visits, function(visit) {
  pre <- visit == "Pre"
  date_candidates <- if (pre) "Pre_CT_Date" else c(paste0("Post_", visit, "_CT_Date"), paste0("Post_", visit, "_CT_Date.1"))
  vfa_candidates <- if (pre) "Pre_L3_VFA_cm2" else paste0("Post_", visit, "_L3_VFA_cm2")
  artifact_candidates <- if (pre) {
    "Pre_L3_Artifact (0=Qualified, (1=Artifact_NonContrast, 2=Artifact_Contrast)"
  } else if (visit == "1st") {
    c("Post_1st_L3_Artifact (0=Qualified, (1=Artifact_NonContrast, 2=Artifact_Contrast)",
      "Post_1st_L2_Artifact (0=Qualified, (1=Artifact_NonContrast, 2=Artifact_Contrast)")
  } else {
    paste0("Post_", visit, "_L3_Artifact (0=Qualified, (1=Artifact_NonContrast, 2=Artifact_Contrast)")
  }
  data.frame(visit, pre, date_col = first_col(nms, date_candidates, FALSE),
             vfa_col = first_col(nms, vfa_candidates, FALSE),
             artifact_col = first_col(nms, artifact_candidates, FALSE))
}))

blocks <- lapply(seq_len(nrow(availability)), function(i) {
  a <- availability[i, ]
  if (any(is.na(a[c("date_col", "vfa_col", "artifact_col")]))) return(NULL)
  data.frame(ID = safe_id(raw[[id_col]]), is_preop = a$pre, visit = a$visit,
             measure_date = parse_date(raw[[a$date_col]]), vfa = safe_num(raw[[a$vfa_col]]),
             artifact = safe_num(raw[[a$artifact_col]]))
})
long_raw <- bind_rows(blocks) %>% inner_join(clinical %>% dplyr::select(ID, surgery_date), by = "ID") %>%
  mutate(time_months = as.numeric(difftime(measure_date, surgery_date, units = "days")) / 30.4375,
         in_window = ifelse(is_preop, time_months >= -3 & time_months <= 0, time_months > 0 & time_months <= 36))
long_qc <- long_raw %>% filter(artifact == 0, in_window) %>% arrange(ID, time_months) %>%
  group_by(ID) %>% mutate(vfa_qc = audit_extreme(vfa)) %>% ungroup()
baseline <- long_qc %>% filter(is_preop, !is.na(vfa_qc)) %>% arrange(ID, desc(time_months)) %>%
  group_by(ID) %>% slice(1) %>% ungroup()
post <- long_qc %>% filter(!is_preop, !is.na(vfa_qc), time_months > 0, time_months <= 24)
post_counts <- post %>% count(ID, name = "n_postoperative_points")
cohort_ids <- baseline %>% dplyr::select(ID) %>% inner_join(post_counts, by = "ID") %>%
  filter(n_postoperative_points >= 2) %>% arrange(ID) %>% mutate(lcmm_subject_id = row_number())
model_data <- bind_rows(baseline, post) %>% filter(ID %in% cohort_ids$ID, time_months <= 24) %>%
  transmute(ID, is_preop, visit, time_months, raw_value = vfa_qc, analysis_value = log1p(vfa_qc)) %>%
  inner_join(cohort_ids, by = "ID")

fixed_formula <- analysis_value ~ splines::ns(time_months, df = 2)
k1 <- do.call(lcmm::hlme, list(
  fixed = fixed_formula, random = ~1, subject = "lcmm_subject_id", ng = 1,
  data = model_data, maxiter = 300, verbose = FALSE
))
k3 <- do.call(lcmm::hlme, list(
  fixed = fixed_formula, mixture = ~splines::ns(time_months, df = 2), random = ~1,
  subject = "lcmm_subject_id", ng = 3, nwg = TRUE, B = k1,
  data = model_data, maxiter = 500, verbose = FALSE
))
if (k3$conv != 1) stop("The K=3 trajectory model did not converge; convergence code: ", k3$conv)

grid <- seq(-3, 24, by = 0.05)
pred0 <- as.data.frame(predictY(k3, newdata = data.frame(time_months = grid),
                               var.time = "time_months", draws = FALSE)$pred)
shape <- bind_rows(lapply(1:3, function(g) data.frame(class = g, time_months = grid, value = pred0[[g]]))) %>%
  group_by(class) %>% summarise(start = value[which.min(abs(time_months))],
                                end = value[which.min(abs(time_months - 24))],
                                mean_level = mean(value), change = end - start, .groups = "drop")
high <- shape$class[which.max(shape$mean_level)]; remaining <- setdiff(shape$class, high)
rebound <- remaining[which.max(shape$change[match(remaining, shape$class)])]
loss <- setdiff(remaining, rebound)
class_map <- setNames(rep(NA_character_, 3), 1:3)
class_map[as.character(high)] <- "High-stable"; class_map[as.character(loss)] <- "Progressive-loss"
class_map[as.character(rebound)] <- "Rebound-increasing"

posterior <- as.data.frame(k3$pprob)
if ("subject" %in% names(posterior)) names(posterior)[names(posterior) == "subject"] <- "lcmm_subject_id"
prob_cols <- grep("^prob", names(posterior), value = TRUE)
prob_matrix <- as.matrix(posterior[, prob_cols, drop = FALSE])
assigned_class <- pmax(1L, pmin(ncol(prob_matrix), safe_num(posterior$class)))
posterior <- posterior %>% mutate(lcmm_subject_id = safe_num(lcmm_subject_id)) %>%
  left_join(cohort_ids, by = "lcmm_subject_id") %>%
  mutate(trajectory = factor(unname(class_map[as.character(class)]), levels = names(colors)),
         assigned_probability = prob_matrix[cbind(seq_len(nrow(prob_matrix)), assigned_class)])

baseline_patient <- baseline %>% filter(ID %in% cohort_ids$ID) %>%
  arrange(ID, desc(time_months)) %>% group_by(ID) %>% slice(1) %>% ungroup() %>%
  transmute(ID, baseline_vfa = vfa_qc, baseline_log_vfa = log1p(vfa_qc))
trajectory_cohort <- posterior %>%
  transmute(ID, trajectory, assigned_probability, n_postoperative_points) %>%
  inner_join(baseline_patient, by = "ID") %>% inner_join(clinical, by = "ID")
landmark <- trajectory_cohort %>% filter(os_months > 24, !is.na(death_event)) %>%
  mutate(os_landmark = os_months - 24,
         age_group = factor(ifelse(age < 65, "<65 years", ">=65 years"), levels = c("<65 years", ">=65 years")),
         stage_group = factor(ifelse(stage %in% c("I", "II"), "Stage I-II", "Stage III"), levels = c("Stage I-II", "Stage III")))

stopifnot(nrow(trajectory_cohort) == 1979, nrow(landmark) == 1844)
class_summary <- trajectory_cohort %>% group_by(trajectory) %>%
  summarise(n = n(), percent = 100 * n() / nrow(trajectory_cohort), AvePP = mean(assigned_probability), .groups = "drop")

fmt_mean <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
fmt_median <- function(x) sprintf("%.1f (%.1f-%.1f)", median(x, na.rm = TRUE),
                                  quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))
fmt_count <- function(data, variable, level) {
  n <- sum(as.character(data[[variable]]) == level, na.rm = TRUE)
  sprintf("%d (%.1f%%)", n, 100 * n / nrow(data))
}
groups <- split(trajectory_cohort, trajectory_cohort$trajectory)
trajectory_characteristics <- data.frame(
  Characteristic = c("Patients", "Average posterior probability", "Age, years", "Sex", "BMI, kg/m2",
                     "Pathologic stage", "", "", "Adjuvant chemotherapy", "", "",
                     "Baseline VFA, cm2", "Postoperative VFA measurements"),
  Level = c("", "", "", "Female", "", "I", "II", "III", "No", "Yes", "Missing", "", ""),
  `High-stable` = c(
    nrow(groups$`High-stable`), sprintf("%.3f", mean(groups$`High-stable`$assigned_probability)),
    fmt_mean(groups$`High-stable`$age), fmt_count(groups$`High-stable`, "sex", "Female"),
    fmt_mean(groups$`High-stable`$bmi), sapply(c("I", "II", "III"), function(z) fmt_count(groups$`High-stable`, "stage", z)),
    sapply(c("No", "Yes"), function(z) fmt_count(groups$`High-stable`, "chemo", z)),
    sprintf("%d (%.1f%%)", sum(is.na(groups$`High-stable`$chemo)), 100 * mean(is.na(groups$`High-stable`$chemo))),
    fmt_median(groups$`High-stable`$baseline_vfa), fmt_median(groups$`High-stable`$n_postoperative_points)
  ),
  `Progressive-loss` = c(
    nrow(groups$`Progressive-loss`), sprintf("%.3f", mean(groups$`Progressive-loss`$assigned_probability)),
    fmt_mean(groups$`Progressive-loss`$age), fmt_count(groups$`Progressive-loss`, "sex", "Female"),
    fmt_mean(groups$`Progressive-loss`$bmi), sapply(c("I", "II", "III"), function(z) fmt_count(groups$`Progressive-loss`, "stage", z)),
    sapply(c("No", "Yes"), function(z) fmt_count(groups$`Progressive-loss`, "chemo", z)),
    sprintf("%d (%.1f%%)", sum(is.na(groups$`Progressive-loss`$chemo)), 100 * mean(is.na(groups$`Progressive-loss`$chemo))),
    fmt_median(groups$`Progressive-loss`$baseline_vfa), fmt_median(groups$`Progressive-loss`$n_postoperative_points)
  ),
  `Rebound-increasing` = c(
    nrow(groups$`Rebound-increasing`), sprintf("%.3f", mean(groups$`Rebound-increasing`$assigned_probability)),
    fmt_mean(groups$`Rebound-increasing`$age), fmt_count(groups$`Rebound-increasing`, "sex", "Female"),
    fmt_mean(groups$`Rebound-increasing`$bmi), sapply(c("I", "II", "III"), function(z) fmt_count(groups$`Rebound-increasing`, "stage", z)),
    sapply(c("No", "Yes"), function(z) fmt_count(groups$`Rebound-increasing`, "chemo", z)),
    sprintf("%d (%.1f%%)", sum(is.na(groups$`Rebound-increasing`$chemo)), 100 * mean(is.na(groups$`Rebound-increasing`$chemo))),
    fmt_median(groups$`Rebound-increasing`$baseline_vfa), fmt_median(groups$`Rebound-increasing`$n_postoperative_points)
  ), check.names = FALSE
)
trajectory_characteristics$`P value` <- c(
  "", "", sprintf("%.3f", kruskal.test(age ~ trajectory, trajectory_cohort)$p.value),
  sprintf("%.3f", chisq.test(table(trajectory_cohort$trajectory, trajectory_cohort$sex))$p.value),
  sprintf("%.3f", kruskal.test(bmi ~ trajectory, trajectory_cohort)$p.value),
  sprintf("%.3f", chisq.test(table(trajectory_cohort$trajectory, trajectory_cohort$stage))$p.value), "", "",
  sprintf("%.3f", chisq.test(table(trajectory_cohort$trajectory, trajectory_cohort$chemo))$p.value), "", "",
  sprintf("%.3f", kruskal.test(baseline_vfa ~ trajectory, trajectory_cohort)$p.value),
  sprintf("%.3f", kruskal.test(n_postoperative_points ~ trajectory, trajectory_cohort)$p.value)
)

# Reset here to reproduce the final confidence bands and imputations.
set.seed(20260814)
pred <- as.data.frame(predictY(k3, newdata = data.frame(time_months = grid),
                              var.time = "time_months", draws = TRUE)$pred)
curves <- bind_rows(lapply(1:3, function(g) data.frame(
  class = g, time_months = grid,
  predicted_log = pred[[paste0("Ypred_class", g)]],
  lower_log = pred[[paste0("lower.Ypred_class", g)]],
  upper_log = pred[[paste0("upper.Ypred_class", g)]]
))) %>% mutate(trajectory = factor(unname(class_map[as.character(class)]), levels = names(colors)),
               predicted_vfa = pmax(0, expm1(predicted_log)), lower_vfa = pmax(0, expm1(lower_log)),
               upper_vfa = pmax(0, expm1(upper_log)))
legend_labels <- setNames(sprintf("%s (n=%d; AvePP=%.2f)", class_summary$trajectory,
                                  class_summary$n, class_summary$AvePP), class_summary$trajectory)
p_traj <- ggplot(curves, aes(time_months, predicted_vfa, colour = trajectory)) +
  geom_ribbon(aes(ymin = lower_vfa, ymax = upper_vfa, fill = trajectory), alpha = 0.16,
              colour = NA, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) + geom_line(linewidth = 1.15) +
  scale_colour_manual(values = colors, labels = legend_labels, drop = FALSE) +
  scale_fill_manual(values = colors, drop = FALSE) +
  scale_x_continuous(breaks = c(-3, 0, 3, 6, 12, 18, 24), limits = c(-3, 24)) +
  labs(x = "Months relative to surgery", y = "Model-estimated VFA (cm2)") + theme_pub(11)
ggsave(file.path(out, "vfa_24m_trajectory.pdf"), p_traj, width = 9, height = 6, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "vfa_24m_trajectory.png"), p_traj, width = 9, height = 6, dpi = 500, bg = "white")

fit_km <- survfit(Surv(os_landmark, death_event) ~ trajectory, data = landmark)
lr <- survdiff(Surv(os_landmark, death_event) ~ trajectory, data = landmark)
lr_p <- pchisq(lr$chisq, df = 2, lower.tail = FALSE)
km <- ggsurvplot(
  fit_km, data = landmark, risk.table = TRUE, risk.table.height = 0.28,
  risk.table.title = "Number at risk", pval = sprintf("Log-rank P = %.3f", lr_p),
  conf.int = FALSE, censor = TRUE, palette = unname(colors), legend.title = NULL,
  legend.labs = sprintf("%s (n=%d)", levels(landmark$trajectory), as.integer(table(landmark$trajectory))),
  xlim = c(0, 96), break.time.by = 12, xlab = "Months after the 24-month landmark",
  ylab = "Overall survival probability", ggtheme = theme_pub(), tables.theme = theme_cleantable(base_size = 9)
)
km$table <- km$table + theme_classic(base_family = "Arial", base_size = 9)
km_grob <- arrangeGrob(km$plot, km$table, heights = c(3.3, 1))
ggsave(file.path(out, "trajectory_km.pdf"), km_grob, width = 9, height = 7.2, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "trajectory_km.png"), km_grob, width = 9, height = 7.2, dpi = 500, bg = "white")

null_cox <- coxph(Surv(os_landmark, death_event) ~ 1, data = landmark)
bh <- basehaz(null_cox, centered = FALSE)
landmark$nelson_aalen <- approx(bh$time, bh$hazard, xout = landmark$os_landmark,
                                method = "constant", f = 0, rule = 2)$y
imp_fit <- glm(chemo_num ~ age + sex + bmi + stage + baseline_log_vfa + trajectory + death_event + nelson_aalen,
               family = binomial(), data = landmark, subset = !is.na(chemo_num))
imp_x <- model.matrix(delete.response(terms(imp_fit)), landmark)
missing_chemo <- which(is.na(landmark$chemo_num))
mi_data <- vector("list", 50)
for (m in 1:50) {
  beta <- MASS::mvrnorm(1, coef(imp_fit), vcov(imp_fit))
  d <- landmark
  d$chemo_num[missing_chemo] <- rbinom(length(missing_chemo), 1,
    plogis(drop(imp_x[missing_chemo, , drop = FALSE] %*% beta)))
  d$chemo <- factor(d$chemo_num, levels = c(0, 1), labels = c("No", "Yes"))
  mi_data[[m]] <- d
}
rng_state_after_mi <- .Random.seed

pool_cox <- function(formula, datasets = mi_data) {
  fits <- lapply(datasets, function(d) coxph(formula, data = d, x = TRUE))
  terms <- names(coef(fits[[1]]))
  table <- bind_rows(lapply(terms, function(term) {
    z <- pool_scalar(vapply(fits, function(f) coef(f)[term], numeric(1)),
                     vapply(fits, function(f) vcov(f)[term, term], numeric(1)))
    data.frame(term, HR = exp(z["estimate"]), CI_low = exp(z["lower"]), CI_high = exp(z["upper"]), p = z["p"])
  }))
  list(fits = fits, table = table)
}
extract_single <- function(fit, term) {
  s <- summary(fit)
  data.frame(term, HR = s$conf.int[term, "exp(coef)"], CI_low = s$conf.int[term, "lower .95"],
             CI_high = s$conf.int[term, "upper .95"], p = s$coefficients[term, "Pr(>|z|)"])
}
terms_traj <- c("trajectoryProgressive-loss", "trajectoryRebound-increasing")
m1 <- coxph(Surv(os_landmark, death_event) ~ trajectory, data = landmark, x = TRUE)
m2 <- pool_cox(Surv(os_landmark, death_event) ~ trajectory + age + sex + bmi + stage + chemo)
m3 <- pool_cox(Surv(os_landmark, death_event) ~ trajectory + age + sex + bmi + stage + chemo + baseline_log_vfa)
cox_results <- bind_rows(
  bind_rows(lapply(terms_traj, function(x) extract_single(m1, x))) %>% mutate(model = "Model 1"),
  m2$table %>% filter(term %in% terms_traj) %>% mutate(model = "Model 2"),
  m3$table %>% filter(term %in% terms_traj) %>% mutate(model = "Model 3")
) %>% mutate(contrast = recode(term,
  `trajectoryProgressive-loss` = "Progressive-loss vs High-stable",
  `trajectoryRebound-increasing` = "Rebound-increasing vs High-stable")) %>%
  dplyr::select(contrast, model, HR, CI_low, CI_high, p)

cc <- landmark %>% filter(!is.na(chemo))
cc_fit <- coxph(Surv(os_landmark, death_event) ~ trajectory + age + sex + bmi + stage + chemo + baseline_log_vfa,
                data = cc, x = TRUE)
cc_results <- bind_rows(lapply(terms_traj, function(x) extract_single(cc_fit, x)))

ph_all <- bind_rows(lapply(seq_along(m3$fits), function(i) {
  z <- cox.zph(m3$fits[[i]], transform = "km", terms = TRUE, global = TRUE)$table
  data.frame(imputation = i, term = rownames(z), chisq = z[, "chisq"], df = z[, "df"], p = z[, "p"], row.names = NULL)
}))
ph_summary <- ph_all %>% group_by(term) %>%
  summarise(n_imputations = n(), n_p_below_0_05 = sum(p < .05), median_p = median(p),
            min_p = min(p), max_p = max(p), .groups = "drop")

write.csv(availability, file.path(out, "ct_column_mapping.csv"), row.names = FALSE)
write.csv(model_data, file.path(out, "trajectory_model_input.csv"), row.names = FALSE)
write.csv(trajectory_cohort, file.path(out, "trajectory_assignments.csv"), row.names = FALSE)
write.csv(class_summary, file.path(out, "trajectory_class_summary.csv"), row.names = FALSE)
write.csv(trajectory_characteristics, file.path(out, "trajectory_characteristics.csv"), row.names = FALSE)
write.csv(data.frame(convergence = k3$conv, log_likelihood = k3$loglik, AIC = k3$AIC,
                     BIC = k3$BIC, patients = nrow(trajectory_cohort), measurements = nrow(model_data)),
          file.path(out, "trajectory_model_summary.csv"), row.names = FALSE)
write.csv(curves, file.path(out, "trajectory_predicted_curves.csv"), row.names = FALSE)
write.csv(cox_results, file.path(out, "trajectory_cox_results.csv"), row.names = FALSE)
write.csv(cc_results, file.path(out, "trajectory_cox_complete_case.csv"), row.names = FALSE)
write.csv(ph_all, file.path(out, "proportional_hazards_tests.csv"), row.names = FALSE)
write.csv(ph_summary, file.path(out, "proportional_hazards_summary.csv"), row.names = FALSE)
saveRDS(k3, file.path(out, "vfa_24m_df2_k3_model.rds"))
saveRDS(list(
  settings = list(trajectory_window = 24, spline_df = 2, classes = 3, landmark = 24,
                  imputations = 50, seed = 20260814),
  trajectory_cohort = trajectory_cohort, landmark = landmark, mi_data = mi_data,
  rng_state_after_mi = rng_state_after_mi, class_summary = class_summary,
  cox_results = cox_results, logrank_p = lr_p
), file.path(out, "analysis_objects.rds"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
print(class_summary)
print(cox_results)
