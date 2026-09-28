options(stringsAsFactors = FALSE, scipen = 999)

# Run after script 02. Both Cox prediction models use identical folds and patients.

output_dir <- Sys.getenv("VFA_OUTPUT_DIR", "results")
local_lib <- Sys.getenv("VFA_R_LIBRARY", "")
if (nzchar(local_lib) && dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
packages <- c("dplyr", "ggplot2", "survival", "MASS", "riskRegression")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(survival); library(MASS); library(riskRegression)
})
if (.Platform$OS.type == "windows") try(windowsFonts(Arial = windowsFont("Arial")), silent = TRUE)

object_file <- file.path(output_dir, "02_trajectory_cox", "analysis_objects.rds")
if (!file.exists(object_file)) stop("Run 02_vfa_trajectory_cox.R first.")
out <- file.path(output_dir, "04_prediction")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

objects <- readRDS(object_file)
pred_data <- objects$landmark %>% arrange(ID)
.Random.seed <- objects$rng_state_after_mi
n_folds <- 10
n_pred_mi <- 20
n_boot <- 2000

make_folds <- function(event, k) {
  out <- integer(length(event))
  for (e in c(0, 1)) {
    idx <- which(event == e)
    out[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  out
}
fold_id <- make_folds(pred_data$death_event, n_folds)
model_names <- c("Clinical + baseline VFA", "Clinical + baseline VFA + trajectory")
oof <- matrix(NA_real_, nrow(pred_data), 2, dimnames = list(NULL, model_names))

imputation_formula <- chemo_num ~ age + sex + bmi + stage + baseline_log_vfa + trajectory
for (fold in seq_len(n_folds)) {
  train_idx <- which(fold_id != fold); validation_idx <- which(fold_id == fold)
  train <- pred_data[train_idx, ]; validation <- pred_data[validation_idx, ]
  imputation_fit <- glm(imputation_formula, family = binomial(), data = train,
                        subset = !is.na(chemo_num))
  terms_no_response <- delete.response(terms(imputation_fit))
  x_train <- model.matrix(terms_no_response, train)
  x_validation <- model.matrix(terms_no_response, validation)
  risk_a <- risk_b <- matrix(NA_real_, nrow(validation), n_pred_mi)

  for (m in seq_len(n_pred_mi)) {
    beta <- MASS::mvrnorm(1, coef(imputation_fit), vcov(imputation_fit))
    completed_train <- train; completed_validation <- validation
    miss_train <- which(is.na(completed_train$chemo_num))
    miss_validation <- which(is.na(completed_validation$chemo_num))
    if (length(miss_train)) {
      completed_train$chemo_num[miss_train] <- rbinom(
        length(miss_train), 1, plogis(drop(x_train[miss_train, , drop = FALSE] %*% beta))
      )
    }
    if (length(miss_validation)) {
      completed_validation$chemo_num[miss_validation] <- rbinom(
        length(miss_validation), 1,
        plogis(drop(x_validation[miss_validation, , drop = FALSE] %*% beta))
      )
    }
    completed_train$chemo <- factor(completed_train$chemo_num, levels = c(0, 1), labels = c("No", "Yes"))
    completed_validation$chemo <- factor(completed_validation$chemo_num, levels = c(0, 1), labels = c("No", "Yes"))

    fit_a <- coxph(Surv(os_landmark, death_event) ~ age + sex + bmi + stage + chemo + baseline_log_vfa,
                    data = completed_train, x = TRUE)
    fit_b <- coxph(Surv(os_landmark, death_event) ~ age + sex + bmi + stage + chemo + baseline_log_vfa + trajectory,
                    data = completed_train, x = TRUE)
    risk_a[, m] <- exp(predict(fit_a, newdata = completed_validation, type = "lp", reference = "zero"))
    risk_b[, m] <- exp(predict(fit_b, newdata = completed_validation, type = "lp", reference = "zero"))
  }
  oof[validation_idx, 1] <- rowMeans(risk_a)
  oof[validation_idx, 2] <- rowMeans(risk_b)
}

harrell_c <- function(index, score) {
  concordance(Surv(pred_data$os_landmark[index], pred_data$death_event[index]) ~ score[index],
              reverse = TRUE)$concordance
}
c_index <- vapply(1:2, function(j) harrell_c(seq_len(nrow(pred_data)), oof[, j]), numeric(1))
bootstrap_index <- replicate(n_boot, sample(seq_len(nrow(pred_data)), replace = TRUE), simplify = FALSE)
bootstrap_c <- sapply(1:2, function(j) vapply(bootstrap_index, harrell_c, numeric(1), score = oof[, j]))

risk_list <- list(
  `Clinical + baseline VFA` = oof[, 1],
  `Clinical + baseline VFA + trajectory` = oof[, 2]
)
auc_score <- riskRegression::Score(
  object = risk_list, formula = Surv(os_landmark, death_event) ~ 1, data = pred_data,
  metrics = "AUC", times = 36, null.model = FALSE, conf.int = .95,
  contrasts = TRUE, cens.method = "ipcw", cens.model = "km", se.fit = TRUE,
  plots = "ROC", verbose = 0
)
auc_table <- as.data.frame(auc_score$AUC$score) %>%
  filter(model != "Null model") %>%
  transmute(Model = as.character(model), AUC = AUC, AUC_low = lower, AUC_high = upper)

performance <- data.frame(
  Model = model_names, N = nrow(pred_data), Events = sum(pred_data$death_event),
  C_index = c_index,
  C_index_low = apply(bootstrap_c, 2, quantile, .025, na.rm = TRUE),
  C_index_high = apply(bootstrap_c, 2, quantile, .975, na.rm = TRUE)
) %>% left_join(auc_table, by = "Model")

delta_c_boot <- bootstrap_c[, 2] - bootstrap_c[, 1]
delta_c <- data.frame(
  comparison = "Trajectory model minus baseline-VFA model",
  delta_C_index = c_index[2] - c_index[1],
  CI_low = quantile(delta_c_boot, .025, na.rm = TRUE),
  CI_high = quantile(delta_c_boot, .975, na.rm = TRUE),
  bootstrap_P = 2 * min(mean(delta_c_boot <= 0), mean(delta_c_boot >= 0))
)
delta_auc <- as.data.frame(auc_score$AUC$contrasts)

write.csv(performance, file.path(out, "prediction_performance.csv"), row.names = FALSE)
write.csv(delta_c, file.path(out, "prediction_delta_c_index.csv"), row.names = FALSE)
write.csv(delta_auc, file.path(out, "prediction_delta_5year_auc.csv"), row.names = FALSE)
write.csv(data.frame(ID = pred_data$ID, fold = fold_id,
                     risk_model_A = oof[, 1], risk_model_B = oof[, 2]),
          file.path(out, "out_of_fold_predictions.csv"), row.names = FALSE)

plot_data <- bind_rows(
  performance %>% transmute(Model, Metric = "Harrell C-index", Estimate = C_index,
                            Lower = C_index_low, Upper = C_index_high),
  performance %>% transmute(Model, Metric = "5-year OS time-dependent AUC", Estimate = AUC,
                            Lower = AUC_low, Upper = AUC_high)
) %>% mutate(Model = factor(Model, levels = rev(model_names)))
p <- ggplot(plot_data, aes(Estimate, Model, colour = Metric)) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = .12, linewidth = .7,
                 position = position_dodge(width = .45)) +
  geom_point(size = 2.5, position = position_dodge(width = .45)) +
  scale_colour_manual(values = c("Harrell C-index" = "#3A6EA5",
                                 "5-year OS time-dependent AUC" = "#C84B31")) +
  scale_x_continuous(limits = c(.5, .8), breaks = seq(.5, .8, .05)) +
  labs(x = "Discrimination estimate", y = NULL) +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = "black"), legend.position = "bottom")
ggsave(file.path(out, "prediction_performance.pdf"), p, width = 9, height = 4.8, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "prediction_performance.png"), p, width = 9, height = 4.8, dpi = 500, bg = "white")

pdf(file.path(out, "five_year_time_dependent_roc.pdf"), width = 6.5, height = 6.5, family = "sans")
plotROC(auc_score, times = 36, legend = TRUE, col = c("#3A6EA5", "#C84B31"), lwd = 2)
dev.off()

capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
print(performance)
print(delta_c)
print(delta_auc)
