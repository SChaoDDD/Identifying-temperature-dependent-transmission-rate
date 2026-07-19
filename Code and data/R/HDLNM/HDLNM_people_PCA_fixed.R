# ============================================================
#        Joint hierarchical DLNM: relative risk analysis
#        Corrected version: PCA coefficients back-transformed
#        to the original cross-basis coefficient space
# ============================================================

# ============================================================
# 1. Load packages
# ============================================================

library(dlnm)
library(splines)
library(glmmTMB)
library(ggplot2)

# ============================================================
# 2. Read data
# ============================================================

# Please modify the following file paths if needed
city_A <- read.csv("C:/Users/XLYC/Desktop/Rprom/HDLNM/IAbeta.csv", header = TRUE)
city_B <- read.csv("C:/Users/XLYC/Desktop/Rprom/HDLNM/ILbeta.csv", header = TRUE)
city_C <- read.csv("C:/Users/XLYC/Desktop/Rprom/HDLNM/INbeta.csv", header = TRUE)

# Standardize column names
colnames(city_A) <- c("time", "beta", "temp")
colnames(city_B) <- c("time", "beta", "temp")
colnames(city_C) <- c("time", "beta", "temp")

# Add city labels
city_A$city <- "IA"
city_B$city <- "IL"
city_C$city <- "IN"

# Combine data
all_data <- rbind(city_A, city_B, city_C)
all_data$city <- factor(all_data$city, levels = c("IA", "IL", "IN"))

# Keep the original row order before merge
all_data$row_id <- seq_len(nrow(all_data))

# Add within-city time index
all_data$time_index <- ave(seq_len(nrow(all_data)), all_data$city, FUN = seq_along)

# ============================================================
# 3. Add population covariate by year
# ============================================================

pop_data <- data.frame(
  city = rep(c("IL", "IA", "IN"), each = 7),
  year = rep(2013:2019, times = 3),
  population = c(
    12895778, 12885092, 12859585, 12821709, 12779893, 12724685, 12667017,
    3093935, 3110643, 3122541, 3133210, 3143734, 3149900, 3159596,
    6570575, 6596019, 6611442, 6637898, 6662068, 6698481, 6731010
  )
)
pop_data$city <- factor(pop_data$city, levels = c("IA", "IL", "IN"))

# Convert time index to date
start_date <- as.Date("2013-09-30")
all_data$date <- start_date + (all_data$time_index - 1)
all_data$year <- as.integer(format(all_data$date, "%Y"))

# Merge population data and restore original order
all_data <- merge(all_data, pop_data, by = c("city", "year"), all.x = TRUE, sort = FALSE)
all_data <- all_data[order(all_data$row_id), ]
rownames(all_data) <- NULL

# Log-transform and standardize population
all_data$log_population <- log(all_data$population)
all_data$log_population_z <- as.numeric(scale(all_data$log_population))

cat("population missing:", sum(is.na(all_data$population)), "\n")
cat("log_population_z missing:", sum(is.na(all_data$log_population_z)), "\n")

# ============================================================
# 4. Data cleaning
# ============================================================

cat("============ Data cleaning ============\n")
cat("Missing value check before cleaning:\n")
cat("  temp missing:", sum(is.na(all_data$temp)), "\n")
cat("  beta missing:", sum(is.na(all_data$beta)), "\n")
cat("  population missing:", sum(is.na(all_data$population)), "\n")
cat("  log_population_z missing:", sum(is.na(all_data$log_population_z)), "\n")
cat("  temp infinite:", sum(is.infinite(all_data$temp)), "\n")
cat("  beta infinite:", sum(is.infinite(all_data$beta)), "\n")
cat("  log_population_z infinite:", sum(is.infinite(all_data$log_population_z)), "\n")

# Remove missing and non-finite rows in original variables
keep_rows <- complete.cases(all_data[, c("temp", "beta", "population", "log_population_z")]) &
  is.finite(all_data$temp) &
  is.finite(all_data$beta) &
  is.finite(all_data$log_population_z)

all_data <- all_data[keep_rows, ]
rownames(all_data) <- NULL

cat("Sample size after cleaning:", nrow(all_data), "\n")

# Gamma model requires a positive response
if (any(all_data$beta <= 0, na.rm = TRUE)) {
  stop("The Gamma(log) model requires beta > 0. Please remove, transform, or use another family for non-positive beta values.")
}

# ============================================================
# 5. Model parameters
# ============================================================

lag_max <- 8
df_var <- 4
df_lag <- 5

temp_median <- median(all_data$temp, na.rm = TRUE)
temp_percentiles <- quantile(all_data$temp, probs = c(0.01, 0.99), na.rm = TRUE)

cat("\n============ Model parameters ============\n")
cat("Maximum lag:", lag_max, "\n")
cat("Reference temperature:", round(temp_median, 2), "C\n")

# ============================================================
# 6. Create cross-basis matrices separately by city
# ============================================================

cat("\n============ Create cross-basis matrix ============\n")

# Use the actual order in the data, not alphabetical factor levels only
city_names <- unique(as.character(all_data$city))
n_cities <- length(city_names)

# Reorder data by city before row-binding city-specific cross-basis matrices.
# This guarantees that cb_combined and all_data are aligned row by row.
all_data_by_city <- do.call(
  rbind,
  lapply(city_names, function(city) all_data[all_data$city == city, ])
)
rownames(all_data_by_city) <- NULL
all_data <- all_data_by_city

cb_list <- list()

for (city in city_names) {
  city_data <- all_data[all_data$city == city, ]
  
  cb_city <- crossbasis(
    x = city_data$temp,
    lag = lag_max,
    argvar = list(fun = "ns", df = df_var, Boundary.knots = temp_percentiles),
    #argvar = list(fun = "ns",knots=equalknots(city_data$temp,fun = "ns", df = 6), Boundary.knots = temp_percentiles),
    arglag = list(fun = "bs", df = df_lag)
  )
  
  cb_list[[city]] <- cb_city
}

cb_combined <- do.call(rbind, cb_list)
n_cb <- ncol(cb_combined)
cb_names <- paste0("cb", seq_len(n_cb))

cat("Cross-basis matrix dimension:", nrow(cb_combined), "x", n_cb, "\n")

# ============================================================
# 7. Remove rows with missing values in the cross-basis matrix
# ============================================================

cat("\n============ Clean cross-basis matrix ============\n")

na_rows <- which(rowSums(is.na(cb_combined)) > 0)
inf_rows <- which(rowSums(!is.finite(cb_combined)) > 0)
problem_rows <- sort(unique(c(na_rows, inf_rows)))

cat("Rows with NA in cross-basis:", length(na_rows), "\n")
cat("Rows with Inf in cross-basis:", length(inf_rows), "\n")
cat("Total rows to remove:", length(problem_rows), "\n")

if (length(problem_rows) > 0) {
  cb_clean <- cb_combined[-problem_rows, , drop = FALSE]
  all_data_clean <- all_data[-problem_rows, , drop = FALSE]
} else {
  cb_clean <- cb_combined
  all_data_clean <- all_data
}

rownames(cb_clean) <- NULL
rownames(all_data_clean) <- NULL

cat("Sample size after cross-basis cleaning:", nrow(all_data_clean), "\n")
cat("Sample size by city:\n")
print(table(all_data_clean$city))

cb_df <- as.data.frame(cb_clean)
colnames(cb_df) <- cb_names
all_data_cb <- cbind(all_data_clean, cb_df)

# ============================================================
# 8. PCA dimension reduction
# ============================================================

cat("\n============ PCA dimension reduction ============\n")

cb_matrix <- as.matrix(cb_df)

cat("Final check - NA count:", sum(is.na(cb_matrix)), "\n")
cat("Final check - non-finite count:", sum(!is.finite(cb_matrix)), "\n")

if (sum(is.na(cb_matrix)) > 0 || sum(!is.finite(cb_matrix)) > 0) {
  stop("The cleaned cross-basis matrix still contains NA or non-finite values.")
}

zero_sd_cols <- which(apply(cb_matrix, 2, sd) == 0)
if (length(zero_sd_cols) > 0) {
  stop(paste0(
    "The following cross-basis columns have zero standard deviation and cannot be scaled in PCA: ",
    paste(colnames(cb_matrix)[zero_sd_cols], collapse = ", ")
  ))
}

pca_result <- prcomp(cb_matrix, center = TRUE, scale. = TRUE)

var_explained <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
n_pc <- which(var_explained >= 0.95)[1]

cat("Selected number of PCs:", n_pc, "\n")
cat("Explained variance:", round(var_explained[n_pc] * 100, 1), "%\n")

pc_scores <- pca_result$x[, seq_len(n_pc), drop = FALSE]
pc_names <- paste0("PC", seq_len(n_pc))
colnames(pc_scores) <- pc_names

all_data_model <- cbind(all_data_cb, as.data.frame(pc_scores))

# ============================================================
# 9. Fit joint hierarchical model
# ============================================================

cat("\n============ Fit joint hierarchical model ============\n")

# Recreate within-city time index after row removal
all_data_model$time_index <- ave(seq_len(nrow(all_data_model)), all_data_model$city, FUN = seq_along)

# Long-term trend control
n_years <- ceiling(max(all_data_model$time_index) / 365)
df_trend <- 12

fixed_part <- paste(pc_names, collapse = " + ")
fixed_formula <- paste(
  "beta ~", fixed_part,
  "+ log_population_z",
  "+ ns(time_index, df =", df_trend, ")"
)

# Random effects: city-specific random intercept plus random slopes for the first PCs
# Use only the first few random slopes to improve convergence.
n_pc_random <- min(3, length(pc_names))
random_part <- paste("(1 +", paste(pc_names[seq_len(n_pc_random)], collapse = " + "), "|| city)")
full_formula <- as.formula(paste(fixed_formula, "+", random_part))

cat("Model formula:\n")
print(full_formula)

cat("\nFitting model...\n")

model_joint <- tryCatch(
  {
    glmmTMB(
      full_formula,
      family = Gamma(link = "log"),
      data = all_data_model,
      control = glmmTMBControl(
        optimizer = optim,
        optArgs = list(method = "BFGS"),
        parallel = 1
      )
    )
  },
  error = function(e) {
    stop(paste("Model fitting failed:", e$message))
  }
)

cat("\nModel fitting completed.\n")
cat("\n============ Model summary ============\n")
print(summary(model_joint))

# ============================================================
# 10. Extract coefficients and back-transform to cross-basis space
# ============================================================

cat("\n============ Coefficient back-transformation ============\n")

fixed_coef <- fixef(model_joint)$cond
pc_coef <- fixed_coef[pc_names]

if (any(is.na(pc_coef))) {
  stop("Some PC coefficients were not found in the fixed effects.")
}

random_coef <- ranef(model_joint)$cond$city
city_names_model <- rownames(random_coef)

# Correct PCA back-transformation:
# prcomp with center=TRUE and scale.=TRUE uses
#   PC = scale(X, center, scale) %*% rotation
# If the fitted PC coefficient vector is gamma, then
#   eta_PC = PC %*% gamma
#          = X %*% [diag(1 / scale) %*% rotation %*% gamma] + constant.
# Therefore, the original cross-basis coefficient vector is
#   beta_cb = diag(1 / scale) %*% rotation[, 1:n_pc] %*% gamma.
# The constant caused by centering is absorbed by the model intercept and is
# irrelevant for relative-risk curves centered by crosspred().
pca_loadings <- pca_result$rotation[, seq_len(n_pc), drop = FALSE]
pca_center <- pca_result$center
pca_scale <- pca_result$scale

if (length(pca_scale) != n_cb) {
  stop("The PCA scale vector length does not match the number of cross-basis columns.")
}
if (any(!is.finite(pca_scale)) || any(pca_scale == 0)) {
  stop("Invalid PCA scale values found. Please check the cross-basis matrix.")
}

# This is the corrected transformation matrix.
# Each row corresponds to one original cross-basis column.
transform_matrix <- sweep(pca_loadings, MARGIN = 1, STATS = pca_scale, FUN = "/")
rownames(transform_matrix) <- cb_names
colnames(transform_matrix) <- pc_names

# Optional constant shift from PCA centering.
# It is not used in crosspred() because RR is computed relative to cen.
pca_intercept_shift <- -sum(pca_center * as.vector(transform_matrix %*% pc_coef))
cat("PCA centering intercept shift, not used in RR prediction:", round(pca_intercept_shift, 6), "\n")

# Pooled cross-basis coefficients
pooled_coef_cb <- as.vector(transform_matrix %*% pc_coef)
names(pooled_coef_cb) <- cb_names

cat("Pooled cross-basis coefficients:\n")
print(round(pooled_coef_cb, 4))

# Fixed-effect covariance transformed to the cross-basis space
vcov_cond <- vcov(model_joint)$cond
vcov_pc <- vcov_cond[pc_names, pc_names, drop = FALSE]
pooled_vcov_cb <- transform_matrix %*% vcov_pc %*% t(transform_matrix)
dimnames(pooled_vcov_cb) <- list(cb_names, cb_names)

# City-specific coefficients.
# Random effects are available only for the first n_pc_random PCs.
# PCs without random slopes receive zero random contribution.
city_coef_list <- list()
city_intercept_shift <- list()

for (i in seq_along(city_names_model)) {
  city_name <- city_names_model[i]
  
  pc_random_effect <- rep(0, length(pc_names))
  names(pc_random_effect) <- pc_names
  
  existing_random_pc <- intersect(pc_names[seq_len(n_pc_random)], colnames(random_coef))
  
  if (length(existing_random_pc) > 0) {
    pc_random_effect[existing_random_pc] <- as.numeric(random_coef[i, existing_random_pc, drop = TRUE])
  }
  
  city_pc_coef <- pc_coef + pc_random_effect
  city_coef_cb <- as.vector(transform_matrix %*% city_pc_coef)
  names(city_coef_cb) <- cb_names
  
  city_coef_list[[city_name]] <- city_coef_cb
  city_intercept_shift[[city_name]] <- -sum(pca_center * city_coef_cb)
}

cat("City-specific PCA centering intercept shifts, not used in RR prediction:\n")
print(round(unlist(city_intercept_shift), 6))

# ============================================================
# 11. Generate predictions
# ============================================================

cat("\n============ Generate predictions ============\n")

temp_pred <- seq(temp_percentiles[1], temp_percentiles[2], length.out = 100)

# Use one fitted cross-basis object as the prediction basis template.
# Its attributes define the same exposure-lag basis used in model construction.
basis_for_prediction <- cb_list[[city_names[1]]]

pred_pooled <- crosspred(
  basis = basis_for_prediction,
  coef = pooled_coef_cb,
  vcov = pooled_vcov_cb,
  model.link = "log",
  at = temp_pred,
  cen = temp_median,
  bylag = 1
)

pred_city <- list()
for (city in city_names_model) {
  pred_city[[city]] <- crosspred(
    basis = basis_for_prediction,
    coef = city_coef_list[[city]],
    vcov = pooled_vcov_cb,
    model.link = "log",
    at = temp_pred,
    cen = temp_median,
    bylag = 1
  )
}

cat("Prediction completed.\n")

# ============================================================
# 12. Base R plots
# ============================================================

cat("\n============ Plot results ============\n")

par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

# Plot 1: pooled cumulative exposure-response curve
plot(
  pred_pooled, "overall",
  xlab = "Temperature (C)",
  ylab = "Cumulative relative risk (RR)",
  main = paste0("Pooled cumulative effect\nReference: ", round(temp_median, 1), " C"),
  col = "red", lwd = 2,
  ci = "area",
  ci.arg = list(col = rgb(1, 0, 0, 0.2))
)
abline(h = 1, lty = 2, col = "gray50")
abline(v = temp_median, lty = 3, col = "blue")

# Plot 2: 3D exposure-lag-response surface
plot(
  pred_pooled, "3d",
  xlab = "Temperature (C)", ylab = "Lag (days)", zlab = "RR",
  main = "Exposure-lag-response surface",
  theta = 40, phi = 30, col = heat.colors(100)
)

# Plot 3: contour plot
plot(
  pred_pooled, "contour",
  xlab = "Temperature (C)", ylab = "Lag (days)",
  main = "Exposure-lag-response contour",
  key.title = title("RR")
)
abline(v = temp_median, lty = 2, col = "white", lwd = 2)

# Plot 4: pooled and city-specific cumulative curves
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

all_rr <- c(pred_pooled$allRRfit, unlist(lapply(pred_city, function(x) x$allRRfit)))
y_range <- range(all_rr, na.rm = TRUE)
y_range[1] <- max(0.3, y_range[1] * 0.9)
y_range[2] <- min(5, y_range[2] * 1.1)

plot(
  pred_pooled$predvar, pred_pooled$allRRfit,
  type = "l",
  xlab = "Temperature (C)",
  ylab = "Cumulative relative risk (RR)",
  main = "Joint hierarchical DLNM: pooled and city-specific effects",
  col = "red", lwd = 3, ylim = y_range
)

polygon(
  c(pred_pooled$predvar, rev(pred_pooled$predvar)),
  c(pred_pooled$allRRlow, rev(pred_pooled$allRRhigh)),
  col = rgb(1, 0, 0, 0.15), border = NA
)

city_colors <- c("#1b9e77", "#d95f02", "#7570b3")
for (i in seq_along(city_names_model)) {
  lines(
    pred_city[[city_names_model[i]]]$predvar,
    pred_city[[city_names_model[i]]]$allRRfit,
    col = city_colors[i], lwd = 2, lty = 2
  )
}

abline(h = 1, lty = 2, col = "gray50")
abline(v = temp_median, lty = 3, col = "gray50")

legend(
  "topright",
  legend = c("Pooled", city_names_model),
  col = c("red", city_colors[seq_along(city_names_model)]),
  lty = c(1, rep(2, length(city_names_model))),
  lwd = c(3, rep(2, length(city_names_model))),
  bty = "n"
)

# ============================================================
# 13. Publication-style cumulative curve, base R version
# ============================================================

city_colors <- c("#1b9e77", "#d95f02", "#7570b3")
temp_seq <- pred_pooled$predvar

all_rr <- c(pred_pooled$allRRfit, unlist(lapply(pred_city, function(x) x$allRRfit)))
y_lim <- c(
  max(0.5, min(all_rr, na.rm = TRUE) * 0.9),
  min(4, max(all_rr, na.rm = TRUE) * 1.1)
)

par(mar = c(5, 4, 4, 2))

plot(
  temp_seq, pred_pooled$allRRfit,
  type = "n",
  xlab = "Temperature (C)",
  ylab = "Cumulative relative risk (RR)",
  main = "Cumulative exposure-response association between temperature and beta",
  ylim = y_lim
)

polygon(
  c(temp_seq, rev(temp_seq)),
  c(pred_pooled$allRRlow, rev(pred_pooled$allRRhigh)),
  col = adjustcolor("red", alpha.f = 0.2), border = NA
)

for (i in seq_along(city_names_model)) {
  lines(
    temp_seq, pred_city[[city_names_model[i]]]$allRRfit,
    col = city_colors[i], lwd = 2, lty = 2
  )
}

lines(temp_seq, pred_pooled$allRRfit, col = "red", lwd = 2.5)

abline(h = 1, lty = 2, col = "gray50")
abline(v = temp_median, lty = 3, col = "gray50")

legend(
  "topright",
  legend = c("Pooled", city_names_model),
  col = c("red", city_colors[seq_along(city_names_model)]),
  lty = c(1, rep(2, length(city_names_model))),
  lwd = 2,
  bty = "n"
)

mtext(paste0("Reference temperature: ", round(temp_median, 1), " C"), side = 3, line = 0, cex = 0.9)

# ============================================================
# 14. Numerical output
# ============================================================

cat("\n============ Key numerical results ============\n")

key_temps <- c(
  quantile(all_data_clean$temp, 0.05),
  temp_median,
  quantile(all_data_clean$temp, 0.95)
)

cat("\nCumulative relative risk, reference temperature:", round(temp_median, 1), "C\n")
cat(rep("-", 60), "\n", sep = "")

for (temp_val in key_temps) {
  idx <- which.min(abs(pred_pooled$predvar - temp_val))
  
  if (abs(temp_val - temp_median) < 0.5) {
    cat(sprintf("\nTemperature %.1f C, reference: RR = 1.00\n", temp_val))
  } else {
    cat(sprintf("\nTemperature %.1f C:\n", temp_val))
    cat(sprintf(
      "  Pooled effect: RR = %.3f (95%% CI: %.3f - %.3f)\n",
      pred_pooled$allRRfit[idx],
      pred_pooled$allRRlow[idx],
      pred_pooled$allRRhigh[idx]
    ))
    for (city in city_names_model) {
      cat(sprintf("  %s: RR = %.3f\n", city, pred_city[[city]]$allRRfit[idx]))
    }
  }
}

# ============================================================
# 15. Export results
# ============================================================

output_results <- data.frame(
  Temperature = pred_pooled$predvar,
  Overall_RR = pred_pooled$allRRfit,
  Overall_RR_low = pred_pooled$allRRlow,
  Overall_RR_high = pred_pooled$allRRhigh
)

for (city in city_names_model) {
  output_results[[paste0(city, "_RR")]] <- pred_city[[city]]$allRRfit
}

cat("\n\nPreview of output_results:\n")
print(head(output_results, 10))

# Uncomment to save the output table
#write.csv(output_results, "C:/Users/XLYC/Desktop/Rprom/HDLNM/sen/totalequal6.csv", row.names = FALSE)

cat("\n\n============ Analysis completed. ============\n")
print(summary(model_joint))
# 三维曲面RR表
rr_mat <- pred_pooled$matRRfit

lag_seq <- seq(pred_pooled$lag[1], pred_pooled$lag[2], by = 1)
colnames(rr_mat) <- paste0("Lag_", lag_seq)

surface_table <- data.frame(
  Temperature = pred_pooled$predvar,
  rr_mat,
  check.names = FALSE
)

#write.csv(surface_table, "C:/Users/XLYC/Desktop/Rprom/HDLNM/DLNM_3D_surface_RR_table.csv", row.names = FALSE)
