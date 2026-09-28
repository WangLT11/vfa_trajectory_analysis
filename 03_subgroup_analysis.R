options(stringsAsFactors = FALSE, scipen = 999)

# Run after scripts 01 and 02. All subgroup models use the saved landmark cohort.

output_dir <- Sys.getenv("VFA_OUTPUT_DIR", "results")
local_lib <- Sys.getenv("VFA_R_LIBRARY", "")
if (nzchar(local_lib) && dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
packages <- c("dplyr", "ggplot2", "survival", "survminer", "gridExtra")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(survival); library(survminer); library(gridExtra)
})
if (.Platform$OS.type == "windows") try(windowsFonts(Arial = windowsFont("Arial")), silent = TRUE)

object_file <- file.path(output_dir, "02_trajectory_cox", "analysis_objects.rds")
cutpoint_file <- file.path(output_dir, "01_preoperative_vfa", "preoperative_vfa_results.csv")
if (!file.exists(object_file)) stop("Run 02_vfa_trajectory_cox.R first.")
if (!file.exists(cutpoint_file)) stop("Run 01_preoperative_vfa_survival.R first.")
out <- file.path(output_dir, "03_subgroup")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

objects <- readRDS(object_file)
landmark <- objects$landmark
mi_data <- objects$mi_data
n_mi <- length(mi_data)
terms_traj <- c("trajectoryProgressive-loss", "trajectoryRebound-increasing")
colors <- c("Progressive-loss" = "#C84B31", "Rebound-increasing" = "#3A6EA5")

pool_scalar <- function(q, u) {
  m <- length(q); qbar <- mean(q); ubar <- mean(u); b <- var(q)
  total <- ubar + (1 + 1 / m) * b
  df <- if (!is.finite(b) || b < 1e-14) Inf else (m - 1) * (1 + ubar / ((1 + 1 / m) * b))^2
  se <- sqrt(total); crit <- if (is.finite(df)) qt(.975, df) else qnorm(.975)
  p <- if (is.finite(df)) 2 * pt(-abs(qbar / se), df) else 2 * pnorm(-abs(qbar / se))
  c(estimate = qbar, lower = qbar - crit * se, upper = qbar + crit * se, p = p)
}
extract_single <- function(fit, term) {
  s <- summary(fit)
  data.frame(HR = s$conf.int[term, "exp(coef)"], CI_low = s$conf.int[term, "lower .95"],
             CI_high = s$conf.int[term, "upper .95"], p = s$coefficients[term, "Pr(>|z|)"])
}
pool_level <- function(group_var, group_level, contrast, covariates) {
  formula <- as.formula(paste("Surv(os_landmark, death_event) ~ trajectory +", paste(covariates, collapse = " + ")))
  q <- u <- n <- events <- numeric(n_mi)
  for (m in seq_len(n_mi)) {
    d <- mi_data[[m]][as.character(mi_data[[m]][[group_var]]) == group_level, , drop = FALSE]
    fit <- coxph(formula, data = d)
    q[m] <- coef(fit)[contrast]; u[m] <- vcov(fit)[contrast, contrast]
    n[m] <- nrow(d); events[m] <- sum(d$death_event)
  }
  z <- pool_scalar(q, u)
  data.frame(N = round(mean(n)), Events = round(mean(events)), HR = exp(z["estimate"]),
             CI_low = exp(z["lower"]), CI_high = exp(z["upper"]), p = z["p"])
}
interaction_p <- function(group_var, contrast, covariates) {
  formula <- as.formula(paste("Surv(os_landmark, death_event) ~ trajectory *", group_var, "+",
                              paste(covariates, collapse = " + ")))
  q <- u <- numeric(n_mi)
  expected <- paste0(contrast, ":", group_var, levels(mi_data[[1]][[group_var]])[2])
  for (m in seq_len(n_mi)) {
    fit <- coxph(formula, data = mi_data[[m]])
    term <- names(coef(fit))[names(coef(fit)) == expected]
    if (!length(term)) term <- names(coef(fit))[grepl(paste0("^", contrast, ":", group_var), names(coef(fit)))]
    q[m] <- coef(fit)[term[1]]; u[m] <- vcov(fit)[term[1], term[1]]
  }
  pool_scalar(q, u)["p"]
}

specs <- list(
  list(label = "Age", var = "age_group", levels = c("<65 years", ">=65 years"),
       covariates = c("age", "sex", "bmi", "stage", "chemo", "baseline_log_vfa")),
  list(label = "Sex", var = "sex", levels = c("Male", "Female"),
       covariates = c("age", "bmi", "stage", "chemo", "baseline_log_vfa")),
  list(label = "Pathologic stage", var = "stage_group", levels = c("Stage I-II", "Stage III"),
       covariates = c("age", "sex", "bmi", "chemo", "baseline_log_vfa"))
)

rows <- list()
for (spec in specs) {
  for (contrast in terms_traj) {
    pint <- interaction_p(spec$var, contrast, spec$covariates)
    for (level in spec$levels) {
      z <- pool_level(spec$var, level, contrast, spec$covariates)
      rows[[length(rows) + 1]] <- cbind(
        data.frame(Contrast = recode(contrast,
          `trajectoryProgressive-loss` = "Progressive-loss",
          `trajectoryRebound-increasing` = "Rebound-increasing"),
          Subgroup = spec$label, Level = level), z, Interaction_P = pint
      )
    }
  }
}
subgroup <- bind_rows(rows)

complete_chemo <- landmark %>% filter(!is.na(chemo))
for (contrast in terms_traj) {
  fit_int <- coxph(Surv(os_landmark, death_event) ~ trajectory * chemo + age + sex + bmi + stage + baseline_log_vfa,
                     data = complete_chemo)
  int_term <- names(coef(fit_int))[grepl(paste0("^", contrast, ":chemo"), names(coef(fit_int)))][1]
  pint <- summary(fit_int)$coefficients[int_term, "Pr(>|z|)"]
  for (level in c("No", "Yes")) {
    d <- complete_chemo %>% filter(as.character(chemo) == level)
    fit <- coxph(Surv(os_landmark, death_event) ~ trajectory + age + sex + bmi + stage + baseline_log_vfa, data = d)
    z <- extract_single(fit, contrast)
    subgroup <- bind_rows(subgroup, data.frame(
      Contrast = recode(contrast,
        `trajectoryProgressive-loss` = "Progressive-loss",
        `trajectoryRebound-increasing` = "Rebound-increasing"),
      Subgroup = "Adjuvant chemotherapy (observed)", Level = level,
      N = nrow(d), Events = sum(d$death_event), z, Interaction_P = pint
    ))
  }
}
write.csv(subgroup, file.path(out, "subgroup_results.csv"), row.names = FALSE)

plot_data <- subgroup %>% filter(Subgroup != "Adjuvant chemotherapy (observed)") %>%
  mutate(label = factor(paste0(Subgroup, ": ", Level), levels = rev(unique(paste0(Subgroup, ": ", Level)))),
         Contrast = factor(Contrast, levels = c("Progressive-loss", "Rebound-increasing")),
         y = as.numeric(label) + ifelse(Contrast == "Progressive-loss", .11, -.11))
p <- ggplot(plot_data, aes(colour = Contrast)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "black", linewidth = .5) +
  geom_segment(aes(x = CI_low, xend = CI_high, y = y, yend = y), linewidth = .65) +
  geom_point(aes(x = HR, y = y), size = 2.3) +
  scale_y_continuous(breaks = seq_along(levels(plot_data$label)), labels = levels(plot_data$label)) +
  scale_x_log10(breaks = c(.5, 1, 2, 4), limits = c(.25, 6)) +
  scale_colour_manual(values = colors) + labs(x = "Hazard ratio (log scale)", y = NULL) +
  theme_classic(base_family = "Arial", base_size = 10.5) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = "black"), legend.position = "bottom")
ggsave(file.path(out, "subgroup_forest_plot.pdf"), p, width = 9, height = 6.5, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "subgroup_forest_plot.png"), p, width = 9, height = 6.5, dpi = 500, bg = "white")

cutpoint <- read.csv(cutpoint_file, check.names = FALSE)$cutpoint_cm2[1]
baseline_strata <- landmark %>%
  mutate(baseline_vfa_group = factor(ifelse(baseline_vfa <= cutpoint, "Low preoperative VFA", "High preoperative VFA"),
                                     levels = c("Low preoperative VFA", "High preoperative VFA")))
baseline_results <- bind_rows(lapply(levels(baseline_strata$baseline_vfa_group), function(level) {
  d <- baseline_strata %>% filter(baseline_vfa_group == level)
  fit <- coxph(Surv(os_landmark, death_event) ~ trajectory, data = d)
  lr <- survdiff(Surv(os_landmark, death_event) ~ trajectory, data = d)
  bind_rows(lapply(terms_traj, function(term) cbind(
    data.frame(stratum = level, N = nrow(d), Events = sum(d$death_event),
               logrank_p = pchisq(lr$chisq, df = 2, lower.tail = FALSE), term = term),
    extract_single(fit, term)
  )))
}))
write.csv(baseline_results, file.path(out, "trajectory_within_preoperative_vfa_strata.csv"), row.names = FALSE)

grDevices::cairo_pdf(file.path(out, "trajectory_km_within_preoperative_vfa_strata.pdf"),
                     width = 8, height = 7, family = "Arial")
for (level in levels(baseline_strata$baseline_vfa_group)) {
  d <- baseline_strata %>% filter(baseline_vfa_group == level)
  fit <- survfit(Surv(os_landmark, death_event) ~ trajectory, data = d)
  lr <- survdiff(Surv(os_landmark, death_event) ~ trajectory, data = d)
  pvalue <- pchisq(lr$chisq, df = 2, lower.tail = FALSE)
  s <- survminer::surv_summary(fit, data = d)
  s$trajectory <- factor(sub(".*=", "", as.character(s$strata)), levels = levels(d$trajectory))
  print(ggplot(s, aes(time, surv, colour = trajectory)) + geom_step(linewidth = .85) +
          annotate("text", x = 3, y = .68, label = sprintf("Log-rank P = %.3f", pvalue),
                   hjust = 0, family = "Arial", size = 3.2) +
          scale_colour_manual(values = c("#2A9D8F", "#C84B31", "#3A6EA5"), drop = FALSE) +
          coord_cartesian(ylim = c(.65, 1.005)) +
          labs(x = "Months after the 24-month landmark", y = "Overall survival probability",
               title = level, colour = NULL) +
          theme_classic(base_family = "Arial", base_size = 10.5) +
          theme(panel.grid = element_blank(), legend.position = "bottom"))
}
dev.off()

capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
print(subgroup)
