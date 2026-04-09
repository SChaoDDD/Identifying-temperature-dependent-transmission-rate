# ============================================================
#        联合层次DLNM模型：流感传染率与温度的相对风险分析
#        （修复版：处理交叉基矩阵中的缺失值）
# ============================================================

# ============================================================
# 第一部分：加载必要的包
# ============================================================

library(dlnm)
library(splines)
library(glmmTMB)
library(ggplot2)

# ============================================================
# 第二部分：读取数据
# ============================================================

# 读取三个城市的数据文件（请修改为实际文件名和路径）
city_A <- read.csv("C:/Users/XLYC/Desktop/流感代码修改/R/IAbeta.csv", header = TRUE)
city_B <- read.csv("C:/Users/XLYC/Desktop/流感代码修改/R/ILbeta.csv", header = TRUE)
city_C <- read.csv("C:/Users/XLYC/Desktop/流感代码修改/R/INbeta.csv", header = TRUE)

# 统一列名
colnames(city_A) <- c("time", "beta", "temp")
colnames(city_B) <- c("time", "beta", "temp")
colnames(city_C) <- c("time", "beta", "temp")

# 添加城市标识
city_A$city <- "City_A"
city_B$city <- "City_B"
city_C$city <- "City_C"

# 合并所有城市数据
all_data <- rbind(city_A, city_B, city_C)
all_data$city <- as.factor(all_data$city)

# 添加时间索引
all_data$time_index <- ave(seq_len(nrow(all_data)), all_data$city, 
                           FUN = seq_along)
# ============================================================
# 第三部分：数据清洗（关键步骤！）
# ============================================================

cat("============ 数据清洗 ============\n")

# 检查原始数据中的缺失值
cat("原始数据缺失值检查:\n")
cat("  temp缺失:", sum(is.na(all_data$temp)), "\n")
cat("  beta缺失:", sum(is.na(all_data$beta)), "\n")

# 检查无穷值
cat("  temp无穷值:", sum(is.infinite(all_data$temp)), "\n")
cat("  beta无穷值:", sum(is.infinite(all_data$beta)), "\n")

# 移除原始数据中的缺失值和无穷值
all_data <- all_data[complete.cases(all_data[, c("temp", "beta")]), ]
all_data <- all_data[is.finite(all_data$temp) & is.finite(all_data$beta), ]

cat("清洗后数据量:", nrow(all_data), "\n")

# ============================================================
# 第四部分：设置模型参数
# ============================================================

lag_max <- 7
df_var <- 3
df_lag <- 4

temp_median <- median(all_data$temp, na.rm = TRUE)
temp_percentiles <- quantile(all_data$temp, probs = c(0.01, 0.99), na.rm = TRUE)

cat("\n============ 模型参数 ============\n")
cat("最大滞后天数:", lag_max, "\n")
cat("基准温度（中位数）:", round(temp_median, 2), "°C\n")

# ============================================================
# 第五部分：为每个城市分别创建交叉基（避免跨城市滞后问题）
# 说明：如果直接对合并数据创建交叉基，会导致城市边界处
#       的滞后计算错误。需要分城市处理。
# ============================================================

cat("\n============ 创建交叉基矩阵 ============\n")

city_names <- levels(all_data$city)
n_cities <- length(city_names)

# 存储各城市的交叉基
cb_list <- list()

for(city in city_names) {
  city_data <- all_data[all_data$city == city, ]
  # 为该城市创建交叉基
  cb_city <- crossbasis(
    x = city_data$temp,
    lag = lag_max,
    argvar = list(fun = "ns", df = df_var, Boundary.knots = temp_percentiles),
    #argvar = list(fun = "ns", df = df_var,knots=equalknots(city_data$temp,fun = "ns", df = 4), Boundary.knots = temp_percentiles),
    arglag = list(fun = "bs", df = df_lag)
  )
  
  cb_list[[city]] <- cb_city
}

# 合并交叉基矩阵
cb_combined <- do.call(rbind, cb_list)
n_cb <- ncol(cb_combined)
cb_names <- paste0("cb", 1:n_cb)

cat("交叉基矩阵维度:", nrow(cb_combined), "×", n_cb, "\n")

# ============================================================
# 第六部分：处理交叉基中的缺失值（关键修复！）
# ============================================================

cat("\n============ 处理交叉基缺失值 ============\n")

# 检查交叉基中的缺失值
na_rows <- which(rowSums(is.na(cb_combined)) > 0)
inf_rows <- which(rowSums(!is.finite(cb_combined)) > 0)
problem_rows <- unique(c(na_rows, inf_rows))

cat("交叉基中含NA的行数:", length(na_rows), "\n")
cat("交叉基中含Inf的行数:", length(inf_rows), "\n")
cat("需要移除的总行数:", length(problem_rows), "\n")

# 移除有问题的行
if(length(problem_rows) > 0) {
  cb_clean <- cb_combined[-problem_rows, ]
  all_data_clean <- all_data[-problem_rows, ]
} else {
  cb_clean <- cb_combined
  all_data_clean <- all_data
}

cat("清洗后数据量:", nrow(all_data_clean), "\n")
cat("各城市数据量:\n")
print(table(all_data_clean$city))

# 将清洗后的交叉基添加到数据框
cb_df <- as.data.frame(cb_clean)
colnames(cb_df) <- cb_names
all_data_cb <- cbind(all_data_clean, cb_df)

# ============================================================
# 第七部分：PCA降维（现在应该没有NA了）
# ============================================================

cat("\n============ PCA降维 ============\n")

# 再次检查确保没有NA
cb_matrix <- as.matrix(cb_df)
cat("最终检查 - NA数量:", sum(is.na(cb_matrix)), "\n")
cat("最终检查 - Inf数量:", sum(!is.finite(cb_matrix)), "\n")

# 进行PCA
pca_result <- prcomp(cb_matrix, center = TRUE, scale. = TRUE)

# 选择解释95%方差的主成分
var_explained <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
n_pc <- which(var_explained >= 0.95)[1]

cat("选择主成分数量:", n_pc, "\n")
cat("解释方差比例:", round(var_explained[n_pc] * 100, 1), "%\n")

# 提取主成分得分
pc_scores <- pca_result$x[, 1:n_pc, drop = FALSE]
pc_names <- paste0("PC", 1:n_pc)
colnames(pc_scores) <- pc_names

# 添加到数据框
all_data_model <- cbind(all_data_cb, as.data.frame(pc_scores))

# ============================================================
# 第八部分：拟合联合层次模型
# ============================================================

cat("\n============ 拟合联合层次模型 ============\n")

# 更新时间索引（因为移除了一些行）
all_data_model$time_index <- ave(seq_len(nrow(all_data_model)), 
                                 all_data_model$city, FUN = seq_along)

# 计算每年的自由度（用于长期趋势控制）
n_years <- ceiling(max(all_data_model$time_index) / 365)
df_trend <-12 # 每年约4个自由度

# 构建模型公式
fixed_part <- paste(pc_names, collapse = " + ")
fixed_formula <- paste("beta ~", fixed_part, "+ ns(time_index, df =", df_trend, ")")

# 随机效应：城市随机截距 + 主成分随机斜率
# || 表示对角协方差结构（提高收敛性）
random_part <- paste("(1 +", paste(pc_names, collapse = " + "), "|| city)")

full_formula <- as.formula(paste(fixed_formula, "+", random_part))

cat("模型公式:\n")
print(full_formula)

# 拟合模型
cat("\n正在拟合模型，请稍候...\n")

model_joint <- tryCatch({
  glmmTMB(
    full_formula,
    family = Gamma(link = "log"),  # Gamma分布适用于正连续数据
    data = all_data_model,
    control = glmmTMBControl(
      optimizer = optim,
      optArgs = list(method = "BFGS"),
      parallel = 1
    )
  )
}, error = function(e) {
  cat("Gamma分布拟合失败，尝试使用高斯分布...\n")
  # 如果Gamma分布失败，尝试对数变换后用高斯分布
  all_data_model$log_beta <- log(all_data_model$beta + 0.001)
  
  formula_gaussian <- as.formula(gsub("beta", "log_beta", 
                                      as.character(full_formula)[2]))
  formula_gaussian <- as.formula(paste("log_beta ~", 
                                       as.character(full_formula)[3]))
  
  glmmTMB(
    formula_gaussian,
    family = gaussian(),
    data = all_data_model,
    control = glmmTMBControl(optimizer = optim, 
                             optArgs = list(method = "BFGS"))
  )
})

cat("\n模型拟合完成！\n")
cat("\n============ 模型摘要 ============\n")
print(summary(model_joint))

# ============================================================
# 第九部分：提取系数并转换回交叉基空间
# ============================================================

cat("\n============ 系数转换 ============\n")

# 提取固定效应系数
fixed_coef <- fixef(model_joint)$cond
pc_coef <- fixed_coef[pc_names]

# 提取随机效应
random_coef <- ranef(model_joint)$cond$city
city_names_model <- rownames(random_coef)

# PCA转换参数
pca_loadings <- pca_result$rotation[, 1:n_pc, drop = FALSE]
pca_scale <- pca_result$scale

# 转换矩阵
transform_matrix <- pca_loadings %*% diag(pca_scale[1:n_pc], nrow = n_pc)

# 总体平均效应
pooled_coef_cb <- as.vector(transform_matrix %*% pc_coef)
names(pooled_coef_cb) <- cb_names

cat("总体平均交叉基系数:\n")
print(round(pooled_coef_cb, 4))

# 协方差矩阵转换
vcov_pc <- vcov(model_joint)$cond[pc_names, pc_names]
pooled_vcov_cb <- transform_matrix %*% vcov_pc %*% t(transform_matrix)

# 城市特异性系数
city_coef_list <- list()
for(i in seq_along(city_names_model)) {
  city_pc_coef <- pc_coef + as.numeric(random_coef[i, pc_names])
  city_coef_list[[city_names_model[i]]] <- as.vector(transform_matrix %*% city_pc_coef)
}

# ============================================================
# 第十部分：生成预测结果
# ============================================================

cat("\n============ 生成预测结果 ============\n")

# 定义预测温度序列
temp_pred <- seq(temp_percentiles[1], temp_percentiles[2], length.out = 100)

# 创建预测用交叉基
cb_pred <- crossbasis(
  x = temp_pred,
  lag = lag_max,
  argvar = list(fun = "ns", df = df_var, Boundary.knots = temp_percentiles),
  #argvar = list(fun = "ns", df = df_var,knots=equalknots(temp_pred,fun = "ns", df = 4), Boundary.knots = temp_percentiles),
  arglag = list(fun = "bs", df = df_lag)
)

# 总体平均效应预测
pred_pooled <- crosspred(
  basis = cb_pred,
  coef = pooled_coef_cb,
  vcov = pooled_vcov_cb,
  model.link = "log",
  at = temp_pred,
  cen = temp_median,
  bylag = 1
)

# 城市特异性预测
pred_city <- list()
for(city in city_names_model) {
  pred_city[[city]] <- crosspred(
    basis = cb_pred,
    coef = city_coef_list[[city]],
    vcov = pooled_vcov_cb,
    model.link = "log",
    at = temp_pred,
    cen = temp_median,
    bylag = 1
  )
}

cat("预测完成！\n")

# ============================================================
# 第十一部分：可视化结果
# ============================================================

cat("\n============ 绑制结果图 ============\n")

# 设置绑图参数
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

# 图1：总体累积暴露-反应曲线
plot(pred_pooled, "overall",
     xlab = "温度 (°C)",
     ylab = "累积相对风险 (RR)",
     main = paste0("总体累积效应曲线\n(基准温度: ", round(temp_median, 1), "°C)"),
     col = "red", lwd = 2,
     ci = "area",
     ci.arg = list(col = rgb(1, 0, 0, 0.2)))
abline(h = 1, lty = 2, col = "gray50")
abline(v = temp_median, lty = 3, col = "blue")

# 图2：3D暴露-滞后-反应曲面
plot(pred_pooled, "3d",
     xlab = "温度 (°C)", ylab = "滞后 (天)", zlab = "RR",
     main = "暴露-滞后-反应曲面",
     theta = 40, phi = 30, col = heat.colors(100))

# 图3：等高线图
plot(pred_pooled, "contour",
     xlab = "温度 (°C)", ylab = "滞后 (天)",
     main = "暴露-滞后-反应等高线图",
     key.title = title("RR"))
abline(v = temp_median, lty = 2, col = "white", lwd = 2)

# ============================================================
# 图5：城市特异性曲线比较
# ============================================================

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

# 确定y轴范围
all_rr <- c(pred_pooled$allRRfit, 
            unlist(lapply(pred_city, function(x) x$allRRfit)))
y_range <- range(all_rr, na.rm = TRUE)
y_range[1] <- max(0.3, y_range[1] * 0.9)
y_range[2] <- min(5, y_range[2] * 1.1)

# 绑制
plot(pred_pooled$predvar, pred_pooled$allRRfit,
     type = "l", xlab = "温度 (°C)", ylab = "累积相对风险 (RR)",
     main = "联合层次模型：总体效应与城市特异性效应",
     col = "red", lwd = 3, ylim = y_range)

polygon(c(pred_pooled$predvar, rev(pred_pooled$predvar)),
        c(pred_pooled$allRRlow, rev(pred_pooled$allRRhigh)),
        col = rgb(1, 0, 0, 0.15), border = NA)

city_colors <- c("#1b9e77", "#d95f02", "#7570b3")
for(i in seq_along(city_names_model)) {
  lines(pred_city[[city_names_model[i]]]$predvar,
        pred_city[[city_names_model[i]]]$allRRfit,
        col = city_colors[i], lwd = 2, lty = 2)
}

abline(h = 1, lty = 2, col = "gray50")
abline(v = temp_median, lty = 3, col = "gray50")

legend("topright",
       legend = c("总体平均效应", city_names_model),
       col = c("red", city_colors[1:length(city_names_model)]),
       lty = c(1, rep(2, length(city_names_model))),
       lwd = c(3, rep(2, length(city_names_model))),
       bty = "n")

# ============================================================
# 第十二部分：ggplot2高质量图
# ============================================================
# ============================================================
#        简洁版绘图代码（基础R，无依赖问题）
# ============================================================

# 颜色设置
city_colors <- c("#1b9e77", "#d95f02", "#7570b3")

# 获取温度序列
temp_seq <- pred_pooled$predvar

# y轴范围
all_rr <- c(pred_pooled$allRRfit, 
            unlist(lapply(pred_city, function(x) x$allRRfit)))
y_lim <- c(max(0.5, min(all_rr, na.rm = TRUE) * 0.9),
           min(4, max(all_rr, na.rm = TRUE) * 1.1))

# 绘图
par(mar = c(5, 4, 4, 2))

plot(temp_seq, pred_pooled$allRRfit,
     type = "n",
     xlab = "温度 (°C)",
     ylab = "累积相对风险 (RR)",
     main = "温度与流感传染率的累积暴露-反应关系",
     ylim = y_lim)

# 总体置信区间
polygon(c(temp_seq, rev(temp_seq)),
        c(pred_pooled$allRRlow, rev(pred_pooled$allRRhigh)),
        col = adjustcolor("red", alpha = 0.2), border = NA)

# 城市曲线
for(i in seq_along(city_names)) {
  lines(temp_seq, pred_city[[city_names[i]]]$allRRfit,
        col = city_colors[i], lwd = 2, lty = 2)
}

# 总体曲线
lines(temp_seq, pred_pooled$allRRfit, col = "red", lwd = 2.5)

# 参考线
abline(h = 1, lty = 2, col = "gray50")
abline(v = temp_median, lty = 3, col = "gray50")

# 图例
legend("topright",
       legend = c("总体效应", city_names),
       col = c("red", city_colors),
       lty = c(1, 2, 2, 2),
       lwd = 2, bty = "n")

# 添加基准温度标注
mtext(paste0("基准温度: ", round(temp_median, 1), "°C"), 
      side = 3, line = 0, cex = 0.9)
# ============================================================
# 第十三部分：数值结果输出
# ============================================================

cat("\n============ 关键数值结果 ============\n")

key_temps <- c(
  quantile(all_data_clean$temp, 0.05),
  temp_median,
  quantile(all_data_clean$temp, 0.95)
)

cat("\n累积相对风险（基准温度:", round(temp_median, 1), "°C）:\n")
cat(rep("-", 60), "\n", sep = "")

for(temp_val in key_temps) {
  idx <- which.min(abs(pred_pooled$predvar - temp_val))
  
  if(abs(temp_val - temp_median) < 0.5) {
    cat(sprintf("\n温度 %.1f°C (基准): RR = 1.00\n", temp_val))
  } else {
    cat(sprintf("\n温度 %.1f°C:\n", temp_val))
    cat(sprintf("  总体效应: RR = %.3f (95%% CI: %.3f - %.3f)\n",
                pred_pooled$allRRfit[idx],
                pred_pooled$allRRlow[idx],
                pred_pooled$allRRhigh[idx]))
    for(city in city_names_model) {
      cat(sprintf("  %s: RR = %.3f\n", city, pred_city[[city]]$allRRfit[idx]))
    }
  }
}


# ============================================================
# 第十四部分：导出结果
# ============================================================

output_results <- data.frame(
  Temperature = pred_pooled$predvar,
  Overall_RR = pred_pooled$allRRfit,
  Overall_RR_low = pred_pooled$allRRlow,
  Overall_RR_high = pred_pooled$allRRhigh
)

for(city in city_names_model) {
  output_results[[paste0(city, "_RR")]] <- pred_city[[city]]$allRRfit
}

cat("\n\n结果数据框预览:\n")
print(head(output_results, 10))

# 保存结果（取消注释以保存）
#write.csv(output_results, "C:/Users/XLYC/Desktop/流感代码修改/R/senany/total_eaualknot4.csv", row.names = FALSE)

cat("\n\n============ 分析完成！ ============\n")
#fixed_effects <- fixef(model_joint)$cond
#random_sd <- sqrt(diag(VarCorr(model_joint)$cond$city))

#relative_importance <- random_sd / abs(fixed_effects[names(random_sd)])