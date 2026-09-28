options(stringsAsFactors = FALSE, scipen = 999)

# Run first. The default input is data/vfa_body_composition.csv.
# Paths can be overridden with VFA_DATA_FILE and VFA_OUTPUT_DIR.

input_file <- Sys.getenv("VFA_DATA_FILE", "data/vfa_body_composition.csv")
input_encoding <- Sys.getenv("VFA_INPUT_ENCODING", "UTF-8")
output_dir <- Sys.getenv("VFA_OUTPUT_DIR", "results")
local_lib <- Sys.getenv("VFA_R_LIBRARY", "")
if (nzchar(local_lib) && dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

packages <- c("data.table", "dplyr", "ggplot2", "survival", "survminer", "gridExtra")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(gridExtra)
})
if (.Platform$OS.type == "windows") try(windowsFonts(Arial = windowsFont("Arial")), silent = TRUE)

out <- file.path(output_dir, "01_preoperative_vfa")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

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
theme_pub <- function(size = 10.5) {
  theme_classic(base_family = "Arial", base_size = size) +
    theme(panel.grid = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black"), legend.position = "bottom",
          plot.background = element_rect(fill = "white", colour = NA))
}

raw <- data.table::fread(input_file, encoding = input_encoding, data.table = FALSE, check.names = FALSE)
nms <- names(raw)
id_col <- first_col(nms, c(" ID", "ID"))
age_col <- first_col(nms, c("Age (years)", "Age"))
surgery_col <- first_col(nms, c("Surgery_Date", "Surgery Date"))
os_col <- first_col(nms, "Overall_Survival_Months")
death_col <- first_col(nms, "Death_Event")
clinical <- data.frame(
  ID = safe_id(raw[[id_col]]), age = safe_num(raw[[age_col]]),
  surgery_date = parse_date(raw[[surgery_col]]), os_months = safe_num(raw[[os_col]]),
  death_raw = safe_num(raw[[death_col]])
) %>%
  mutate(
    death_event = case_when(death_raw == 0 ~ 0L, death_raw == 1 ~ 1L, TRUE ~ NA_integer_)
  ) %>%
  filter(!is.na(ID), !is.na(age), age < 85, !is.na(surgery_date))

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
  data.frame(ID = safe_id(raw[[id_col]]), is_preop = a$pre,
             measure_date = parse_date(raw[[a$date_col]]), vfa = safe_num(raw[[a$vfa_col]]),
             artifact = safe_num(raw[[a$artifact_col]]))
})
long_qc <- bind_rows(blocks) %>% inner_join(clinical, by = "ID") %>%
  mutate(time_months = as.numeric(difftime(measure_date, surgery_date, units = "days")) / 30.4375,
         in_window = ifelse(is_preop, time_months >= -3 & time_months <= 0,
                            time_months > 0 & time_months <= 36)) %>%
  filter(artifact == 0, in_window) %>% arrange(ID, time_months) %>%
  group_by(ID) %>% mutate(vfa_qc = audit_extreme(vfa)) %>% ungroup() %>%
  filter(is_preop, !is.na(vfa_qc)) %>% arrange(ID, desc(time_months)) %>%
  group_by(ID) %>% slice(1) %>% ungroup() %>%
  filter(!is.na(os_months), os_months > 0, !is.na(death_event)) %>%
  transmute(ID, time = os_months, event = death_event, baseline_vfa = vfa_qc)
patients <- long_qc

cut_object <- surv_cutpoint(
  patients, time = "time", event = "event", variables = "baseline_vfa",
  minprop = 0.10, progressbar = FALSE
)
cutpoint <- as.numeric(cut_object$cutpoint$cutpoint)
categorized <- surv_categorize(cut_object)
class(categorized) <- "data.frame"
categorized <- categorized %>%
  mutate(group = factor(tolower(baseline_vfa), levels = c("low", "high"), labels = c("Low", "High")))
fit <- survfit(Surv(time, event) ~ group, data = categorized)
cox <- coxph(Surv(time, event) ~ group, data = categorized)
lr <- survdiff(Surv(time, event) ~ group, data = categorized)
logrank_p <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)

cut_stats <- data.frame(cut = cut_object$baseline_vfa$cuts, statistic = cut_object$baseline_vfa$stats)
p_cut <- ggplot(cut_stats, aes(cut, statistic)) +
  geom_line(data = filter(cut_stats, cut <= cutpoint), colour = "#2166AC", linewidth = 0.9) +
  geom_line(data = filter(cut_stats, cut >= cutpoint), colour = "#D73027", linewidth = 0.9) +
  geom_vline(xintercept = cutpoint, linetype = "dashed", linewidth = 0.55) +
  annotate("text", x = cutpoint, y = max(cut_stats$statistic),
           label = sprintf("Cutpoint = %.2f cm2", cutpoint), hjust = -0.05, vjust = 1.1,
           family = "Arial", size = 3.2) +
  labs(x = "Preoperative VFA (cm2)", y = "Standardized log-rank statistic",
       title = "A  Cutpoint diagnostic") + theme_pub()

km <- ggsurvplot(
  fit, data = categorized, risk.table = TRUE, risk.table.height = 0.27,
  risk.table.title = "Number at risk", pval = TRUE, conf.int = FALSE,
  censor = TRUE, palette = c("#2166AC", "#D73027"), legend.title = NULL,
  legend.labs = c(sprintf("Low VFA (<= %.2f cm2)", cutpoint),
                  sprintf("High VFA (> %.2f cm2)", cutpoint)),
  xlim = c(0, 120), break.time.by = 24, xlab = "Time from surgery (months)",
  ylab = "Overall survival probability", title = "B  Overall survival by preoperative VFA",
  ggtheme = theme_pub(), tables.theme = theme_cleantable(base_size = 9)
)
km$table <- km$table + theme_classic(base_family = "Arial", base_size = 9)
combined <- arrangeGrob(p_cut, arrangeGrob(km$plot, km$table, heights = c(3.2, 1)), ncol = 2,
                        widths = c(0.9, 1.35))
ggsave(file.path(out, "preoperative_vfa_cutpoint_and_km.pdf"), combined,
       width = 13, height = 6.5, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "preoperative_vfa_cutpoint_and_km.png"), combined,
       width = 13, height = 6.5, dpi = 500, bg = "white")

ci <- summary(cox)$conf.int[1, ]
result <- data.frame(
  n = nrow(categorized), events = sum(categorized$event), cutpoint_cm2 = cutpoint,
  low_n = sum(categorized$group == "Low"), low_events = sum(categorized$event[categorized$group == "Low"]),
  high_n = sum(categorized$group == "High"), high_events = sum(categorized$event[categorized$group == "High"]),
  high_vs_low_HR = ci["exp(coef)"], CI_low = ci["lower .95"], CI_high = ci["upper .95"],
  cox_p = summary(cox)$coefficients[1, "Pr(>|z|)"], logrank_p = logrank_p
)
write.csv(result, file.path(out, "preoperative_vfa_results.csv"), row.names = FALSE)
write.csv(patients, file.path(out, "preoperative_vfa_patient_data.csv"), row.names = FALSE)
write.csv(categorized, file.path(out, "preoperative_vfa_analysis_data.csv"), row.names = FALSE)
saveRDS(cut_object, file.path(out, "preoperative_vfa_cutpoint_model.rds"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
print(result)
