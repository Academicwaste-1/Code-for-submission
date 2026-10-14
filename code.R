
rm(list = ls()) # 清空环境
# -------------------- 0. 只需要修改这两个路径 --------------------
INPUT_FILE <- "C:/Users/zcl/Desktop/data.xlsx"
OUTPUT_DIR <- "C:/Users/zcl/Desktop/Reviewer10_Results_2Metrics"

# 已插补表格中的 BAA 和分组默认沿用，不自动重算，确保与原论文一致。
# 仅使用分析表中已有的 BAA 与分组，避免改变原论文结果。
USE_SAVED_BAA <- TRUE

# MICE 需要未插补的原始分析表，默认关闭；原插补参数见下方函数。
RUN_MICE <- FALSE
PRE_IMPUTATION_FILE <- ""  # 只有明确开启 RUN_MICE 后才填写。

# 是否计算 RCS 图（使用 stats::glm 构建三节点样条模型）。
RUN_RCS <- TRUE


# -------------------- 1. 环境与辅助函数 --------------------
required_packages <- c("readxl", "broom", "ggplot2")
if (RUN_MICE) required_packages <- c(required_packages, "mice")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("缺少 R 包：", paste(missing_packages, collapse = ", "),
       "\n请先运行：install.packages(c(",
       paste(sprintf('"%s"', missing_packages), collapse = ", "), "))")
}
if (!USE_SAVED_BAA) {
  stop("为避免在未核实 CRP 单位前覆盖原结果，本整合版只使用已保存的 BAA；请设 USE_SAVED_BAA <- TRUE。")
}

if (!file.exists(INPUT_FILE)) {
  stop("找不到文件：", INPUT_FILE,
       "\n请确认 Excel 存在，或修改脚本开头的 INPUT_FILE。")
}
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

read_analysis_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    df <- as.data.frame(readxl::read_excel(path, na = c("", "NA"), sheet = 1))
  } else if (ext == "csv") {
    df <- read.csv(path, na.strings = c("", "NA"), check.names = FALSE,
                   fileEncoding = "UTF-8")
  } else stop("仅支持 xlsx、xls 或 csv 数据文件。")
  # 作者的原 Excel 表头含：BioAgeAccel_2011sifenweishu[末尾空格]
  names(df) <- make.unique(trimws(names(df)))
  df
}

write_out <- function(df, name) {
  dest <- file.path(OUTPUT_DIR, name)
  utils::write.csv(df, dest, row.names = FALSE, na = "")
  message("已生成：", dest)
  invisible(dest)
}

check_columns <- function(df, vars) {
  miss <- setdiff(vars, names(df))
  if (length(miss)) {
    stop("缺少这些分析变量：", paste(miss, collapse = ", "),
         "\n请检查 Excel 表头是否与原始分析文件一致。")
  }
}

as_numeric_safe <- function(x) suppressWarnings(as.numeric(as.character(x)))

base_covars <- c("nianling2011", "hukou", "wenhua", "hunyin", "xingbie",
                 "zibaogaoxueya2011", "zibaotangniaobing2011",
                 "xiyan2011", "yinjiu2011", "shousuoya", "shuzhangya", "BMI")
cat_covars <- c("hukou", "wenhua", "hunyin", "xingbie", "zibaogaoxueya2011",
                "zibaotangniaobing2011", "xiyan2011", "yinjiu2011")
yvar <- "shifouxuezhiyichang2015"

# -------------------- 2. 读取并检查 --------------------
message("\n[1/6] 正在读取分析数据 ...")
dat <- read_analysis_file(INPUT_FILE)
check_columns(dat, c(yvar, base_covars, "nianling2015", "BioAge_2011",
                     "BioAge_2015", "BioAgeAccel_2011", "BioAgeAccel_2015",
                     "BAA_Transition", "BAA_Slope", "BAA_Cumulative"))
for (v in c(yvar, "nianling2011", "nianling2015", "shousuoya", "shuzhangya",
            "BMI", "BioAge_2011", "BioAge_2015", "BioAgeAccel_2011",
            "BioAgeAccel_2015", "BAA_Slope", "BAA_Cumulative")) {
  dat[[v]] <- as_numeric_safe(dat[[v]])
}
if (!all(stats::na.omit(unique(dat[[yvar]])) %in% 0:1)) {
  stop("结局变量不是 0/1，请检查 shifouxuezhiyichang2015。")
}
for (v in cat_covars) dat[[v]] <- factor(dat[[v]])

# 优先沿用原始已保存分组；若四分位数列不存在，才用连续数据新建。
quartile_columns <- c("BioAgeAccel_2011sifenweishu", "BAA_Slopesifenweishu",
                      "BAA_Cumulativesifenweishu")
quartile_sources <- c("BioAgeAccel_2011", "BAA_Slope", "BAA_Cumulative")
quartile_fallback <- function(x) {
  if (all(is.na(x))) stop("连续暴露全部缺失，无法划分四分位数。")
  br <- as.numeric(stats::quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE))
  if (length(unique(br)) != 5L) stop("四分位数存在重复分位点，请手动核实。")
  as.character(cut(x, breaks = br, include.lowest = TRUE, labels = as.character(1:4)))
}
for (i in seq_along(quartile_columns)) {
  v <- quartile_columns[i]
  if (!v %in% names(dat)) {
    warning("原表缺少分组列 ", v, "；按连续变量重新计算，结果可能有微小差异。")
    dat[[v]] <- quartile_fallback(dat[[quartile_sources[i]]])
  }
  dat[[v]] <- factor(trimws(as.character(dat[[v]])), levels = as.character(1:4))
}
dat$BAA_Transition <- factor(trimws(as.character(dat$BAA_Transition)),
                             levels = as.character(1:4))
dat$Baseline_BAA_status <- factor(
  ifelse(is.na(dat$BioAgeAccel_2011), NA_character_,
         ifelse(dat$BioAgeAccel_2011 > 0, "2", "1")), levels = c("1", "2")
)

qc <- data.frame(
  Item = c("Total participants", "Incident dyslipidemia", "Outcome missing",
           "2011 accelerated (>0)", "2015 accelerated (>0)",
           "Original transition code 1", "Original transition code 2",
           "Original transition code 3", "Original transition code 4"),
  Value = c(nrow(dat), sum(dat[[yvar]] == 1, na.rm = TRUE), sum(is.na(dat[[yvar]])),
            sum(dat$BioAgeAccel_2011 > 0, na.rm = TRUE),
            sum(dat$BioAgeAccel_2015 > 0, na.rm = TRUE),
            vapply(as.character(1:4), function(k) sum(dat$BAA_Transition == k, na.rm = TRUE),
                   integer(1)))
)
write_out(qc, "01_Data_QC.csv")
print(qc, row.names = FALSE)

# -------------------- 3. LightBioAge / BAA 计算函数 --------------------
# 保留原始方法，默认不覆盖数据中已经保存的变量。

# 函数独立保留，以便公开代码时展示公式、年龄残差化和纵向指标算法。
# 默认不覆盖从 Excel 读取的变量，保持与原稿分析版本一致。
rebuild_lightbioage_BAA <- function(df, crp_unit = "mg_dL", interval = 4) {
  stopifnot(crp_unit %in% c("mg_dL", "mg_L"))
  vars <- c("nianling2011", "nianling2015", "newcrea2011", "newglu2011",
            "newcrp2011", "bl_crea2015", "bl_glu2015", "bl_crp2015")
  check_columns(df, vars)
  for (v in vars) df[[v]] <- as_numeric_safe(df[[v]])
  for (v in c("newcrp2011", "bl_crp2015")) {
    if (any(df[[v]] <= 0, na.rm = TRUE))
      stop("存在非正数 CRP，请先核查检出限处理，不应直接取 log：", v)
  }
  calc_age <- function(age, cr, gluc, crp) {
    crp <- if (crp_unit == "mg_dL") crp / 10 else crp
    age + 8.3313 * cr + 0.8270 * (gluc / 18) +
      5.7305 * log(crp) - 13.5298
  }
  df$BioAge_2011 <- calc_age(df$nianling2011, df$newcrea2011,
                             df$newglu2011, df$newcrp2011)
  df$BioAge_2015 <- calc_age(df$nianling2015, df$bl_crea2015,
                             df$bl_glu2015, df$bl_crp2015)
  df$BioAgeAccel_2011 <- stats::residuals(stats::lm(
    BioAge_2011 ~ nianling2011, data = df, na.action = stats::na.exclude))
  df$BioAgeAccel_2015 <- stats::residuals(stats::lm(
    BioAge_2015 ~ nianling2015, data = df, na.action = stats::na.exclude))

  df$BAA_Slope <- (df$BioAgeAccel_2015 - df$BioAgeAccel_2011) / interval
  df$BAA_Cumulative <- (df$BioAgeAccel_2015 + df$BioAgeAccel_2011) / 2 * interval

  status_2011 <- ifelse(is.na(df$BioAgeAccel_2011), NA_integer_,
                        ifelse(df$BioAgeAccel_2011 > 0, 2L, 1L))
  status_2015 <- ifelse(is.na(df$BioAgeAccel_2015), NA_integer_,
                        ifelse(df$BioAgeAccel_2015 > 0, 2L, 1L))
  df$BAA_Transition <- ifelse(is.na(status_2011) | is.na(status_2015), NA_integer_,
                              ifelse(status_2011 == 1L & status_2015 == 1L, 1L,
                              ifelse(status_2011 == 1L & status_2015 == 2L, 2L,
                              ifelse(status_2011 == 2L & status_2015 == 1L, 3L, 4L))))
  df$BioAgeAccel_2011sifenweishu <- quartile_fallback(df$BioAgeAccel_2011)
  df$BAA_Slopesifenweishu <- quartile_fallback(df$BAA_Slope)
  df$BAA_Cumulativesifenweishu <- quartile_fallback(df$BAA_Cumulative)
  df
}
# -------------------- 4. 多变量 Logistic 回归 --------------------
message("\n[2/6] 正在运行多变量 Logistic 回归 ...")
model_results <- function(df, exposure, ref = "1", adjust_baseline_BAA = FALSE) {
  terms <- base_covars
  if (adjust_baseline_BAA) terms <- c(terms, "BioAgeAccel_2011")
  check_columns(df, c(exposure, yvar, terms))
  d <- df
  if (is.factor(d[[exposure]])) {
    d[[exposure]] <- stats::relevel(d[[exposure]], ref = ref)
  }
  fm <- stats::reformulate(c(exposure, terms), response = yvar)
  fit <- stats::glm(fm, data = d, family = stats::binomial())
  # 与先前脚本一致：broom::tidy(conf.int=TRUE, exponentiate=TRUE)
  result <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  result <- result[startsWith(result$term, exposure),
                   c("term", "estimate", "conf.low", "conf.high", "p.value"),
                   drop = FALSE]
  if (!nrow(result)) stop("模型中未找到暴露系数：", exposure)
  mf <- stats::model.frame(fit)
  result$Exposure <- exposure
  result$Reference <- if (is.factor(d[[exposure]])) ref else NA_character_
  result$N_model <- stats::nobs(fit)
  result$N_cases <- sum(stats::model.response(mf) == 1, na.rm = TRUE)
  result$Adjustment <- if (adjust_baseline_BAA)
    "Main covariates + baseline BAA" else "Main covariates"
  result
}

model_specs <- list(
  list(x = "Baseline_BAA_status", ref = "1", base = FALSE),
  list(x = "BioAgeAccel_2011sifenweishu", ref = "1", base = FALSE),
  list(x = "BAA_Transition", ref = "1", base = FALSE),
  list(x = "BAA_Slopesifenweishu", ref = "1", base = TRUE),
  list(x = "BAA_Cumulativesifenweishu", ref = "1", base = TRUE)
)
results_list <- lapply(model_specs, function(s) {
  message("  模型：", s$x)
  model_results(dat, s$x, ref = s$ref, adjust_baseline_BAA = s$base)
})
logistic_results <- do.call(rbind, results_list)
rownames(logistic_results) <- NULL
write_out(logistic_results, "03_Logistic_Results.csv")
print(logistic_results, row.names = FALSE)

message("  另计算 BAA transition 以第 4 组为参照 ...")
reverse_transition <- model_results(dat, "BAA_Transition", ref = "4",
                                    adjust_baseline_BAA = FALSE)
write_out(reverse_transition, "04_Transition_Ref4.csv")
print(reverse_transition, row.names = FALSE)


message("\n[3/6] 正在进行三节点 RCS 分析 ...")

rcs_nonlinear <- function(x, knots) {
  stopifnot(length(knots) == 3L, all(is.finite(knots)),
            knots[1] < knots[2], knots[2] < knots[3])
  a <- knots[1]; b <- knots[2]; c <- knots[3]
  h <- function(z) pmax(z, 0)^3
  (h(x-a) - h(x-b)*(c-a)/(c-b) + h(x-c)*(b-a)/(c-b))/(c-a)^2
}

run_single_rcs <- function(df, exposure, file_out) {
  covars <- c(base_covars, "BioAgeAccel_2011")
  usevars <- unique(c(yvar, exposure, covars))
  check_columns(df, usevars)
  d <- df[, usevars, drop = FALSE]
  d[[exposure]] <- as_numeric_safe(d[[exposure]])
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nrow(d) < 30L) stop("完整记录少于 30：", exposure)
  if (length(unique(d[[exposure]])) < 5L) stop("连续暴露缺乏变异：", exposure)
  if (length(unique(d[[yvar]])) != 2L) stop("结局不包含两个类别：", exposure)

  x <- d[[exposure]]
  knots <- as.numeric(stats::quantile(x, probs = c(.10, .50, .90), na.rm = TRUE))
  if (length(unique(knots)) < 3L) stop("分位数节点重复，无法拟合：", exposure)
  median_ref <- stats::median(x)
  d$.rcs_x <- x
  d$.rcs_nl <- rcs_nonlinear(x, knots)
  form <- stats::reformulate(c(".rcs_x", ".rcs_nl", covars), response = yvar)
  fit <- stats::glm(form, family = stats::binomial(), data = d,
                    control = stats::glm.control(maxit = 100))
  if (!fit$converged) stop("RCS glm 未收敛：", exposure)

  name2 <- c(".rcs_x", ".rcs_nl")
  beta <- stats::coef(fit)[name2]
  vc <- stats::vcov(fit)[name2, name2, drop = FALSE]
  if (any(!is.finite(beta)) || any(!is.finite(vc)))
    stop("RCS 回归系数/协方差存在非有限值：", exposure)
  # Wald test: H0 beta_linear=beta_nonlinear=0, 2 df.
  wald_overall <- as.numeric(t(beta) %*% solve(vc, beta))
  p_overall <- stats::pchisq(wald_overall, df = 2, lower.tail = FALSE)
  # Wald test: H0 beta_nonlinear=0, 1 df.
  wald_nonlinear <- beta[2]^2 / vc[2,2]
  p_nonlinear <- stats::pchisq(wald_nonlinear, df = 1, lower.tail = FALSE)

  # 预测只画第 1--99 百分位，避免个别极端观测使图像横轴不可读；
  # 回归模型本身仍使用全部可分析样本。
  draw_range <- as.numeric(stats::quantile(x, probs = c(.01, .99)))
  xgrid <- sort(unique(c(seq(draw_range[1], draw_range[2], length.out = 500),
                         median_ref)))
  dx <- xgrid - median_ref
  dz <- rcs_nonlinear(xgrid, knots) - rcs_nonlinear(median_ref, knots)
  log_or <- beta[1]*dx + beta[2]*dz
  var_log_or <- vc[1,1]*dx^2 + 2*vc[1,2]*dx*dz + vc[2,2]*dz^2
  se <- sqrt(pmax(var_log_or, 0))
  pred <- data.frame(x = xgrid, OR = exp(log_or),
                     Lower = exp(log_or - 1.96*se),
                     Upper = exp(log_or + 1.96*se))
  if (any(!is.finite(as.matrix(pred))))
    stop("预测 OR/置信区间非有限，检查极端值：", exposure)

  axis_labels <- c(
    BAA_Slope = "BAA slope (years/year)",
    BAA_Cumulative = "Cumulative BAA (year-years)"
  )
  label <- unname(axis_labels[exposure])
  if (is.na(label)) label <- exposure
  p <- ggplot2::ggplot(pred, ggplot2::aes(x = x, y = OR)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper),
                         fill = "grey75", alpha = 0.35) +
    ggplot2::geom_line(color = "#C34343", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
    ggplot2::geom_point(data = data.frame(x = median_ref, OR = 1),
                        inherit.aes = FALSE,
                        ggplot2::aes(x = x, y = OR), color = "#C34343", size = 2.4) +
    ggplot2::annotate("text", x = median_ref, y = 1,
                      label = paste0("Reference = ", round(median_ref, 2)),
                      vjust = 2, size = 3.2) +
    ggplot2::labs(x = label, y = "Adjusted OR (95% CI)",
                  title = paste("Restricted cubic spline:", label),
                  subtitle = sprintf("Wald P overall = %.3g; P nonlinear = %.3g",
                                     p_overall, p_nonlinear)) +
    ggplot2::theme_classic(base_size = 12)
  ggplot2::ggsave(filename = file_out, plot = p, width = 7, height = 5,
                  units = "in", dpi = 400)
  message("    RCS 成功：", exposure, "; N=", nrow(d),
          "; P非线性=", signif(p_nonlinear, 5),
          "; 图像=", file_out)

  data.frame(Exposure = exposure, N = nrow(d),
             Reference_median = median_ref,
             Knot_10th = knots[1], Knot_50th = knots[2], Knot_90th = knots[3],
             P_overall_Wald = p_overall, P_nonlinear_Wald = p_nonlinear,
             Test = "Wald chi-square; 3-knot restricted cubic spline")
}

if (RUN_RCS) {
  targets <- c("BAA_Slope", "BAA_Cumulative")
  # 不再把实际错误转换为 NA：失败会明确指出是哪一项分析出错。
  rcs_rows <- lapply(targets, function(x) {
    message("  处理暴露指标：", x)
    run_single_rcs(dat, x, file.path(OUTPUT_DIR, paste0("RCS_", x, ".png")))
  })
  rcs_results <- do.call(rbind, rcs_rows)
  rownames(rcs_results) <- NULL
  write_out(rcs_results, "05_RCS_Results.csv")
  print(rcs_results, row.names = FALSE)
  # 另外保存成 RDS，方便 RStudio 随时查看。
  saveRDS(rcs_results, file.path(OUTPUT_DIR, "05_RCS_Results.rds"))
} else {
  message("  RUN_RCS=FALSE，本次未运行 RCS。")
}

# -------------------- 6. E-value（忠实展示原脚本的计算方式） --------------------
message("\n[4/6] 正在计算原脚本 E-value（仅用于复核） ...")
# 原始脚本直接将 OR 带入 RR 的 E-value 公式；对非罕见结局需另行审查。
calc_original_E <- function(x) {
  if (is.na(x) || !is.finite(x) || x <= 0) return(NA_real_)
  rr <- if (x >= 1) x else 1 / x
  rr + sqrt(rr * (rr - 1))
}
evalue_results <- logistic_results

evalue_results$Original_E_from_OR <- vapply(
  evalue_results$estimate, calc_original_E, numeric(1))
evalue_results$Original_E_low <- vapply(
  evalue_results$conf.low, calc_original_E, numeric(1))
evalue_results$Original_E_high <- vapply(
  evalue_results$conf.high, calc_original_E, numeric(1))
evalue_results$E_CI_nearest_null <- mapply(function(est, lcl, ucl) {
  if (anyNA(c(est, lcl, ucl))) return(NA_real_)
  if (lcl <= 1 && ucl >= 1) return(1)
  if (est >= 1) calc_original_E(lcl) else calc_original_E(ucl)
}, evalue_results$estimate, evalue_results$conf.low, evalue_results$conf.high)
write_out(evalue_results, "06_Evalue_Original_Formula_CHECK.csv")
print(evalue_results[, c("Exposure", "term", "Original_E_from_OR",
                         "E_CI_nearest_null")], row.names = FALSE)

# -------------------- 7. R 环境与版本信息 --------------------
message("\n[5/6] 正在记录 R 环境与软件包版本 ...")
writeLines(capture.output(utils::sessionInfo()),
           file.path(OUTPUT_DIR, "07_sessionInfo.txt"))
packages_to_report <- c("mice", "rms", "ggplot2", "readxl", "broom",
                        "dplyr", "lmtest", "EValue")
package_versions <- data.frame(
  Package = packages_to_report,
  Installed_Version = vapply(packages_to_report, function(p) {
    if (requireNamespace(p, quietly = TRUE))
      as.character(utils::packageVersion(p)) else NA_character_
  }, character(1)), row.names = NULL
)
write_out(package_versions, "08_Package_Versions.csv")

# -------------------- 8. MICE 原脚本：只有明确请求时才运行 --------------------
run_original_mice <- function(pre_imputation_df) {
  if (!requireNamespace("mice", quietly = TRUE)) stop("请安装 mice 包。")
  # 来源：mice(data, m=5, maxit=50, seed=826, method='rf')。
  # 仅保留原流程；完整多重插补统计推断需分别拟合 5 套并按 Rubin 法则 pool。
  mice::mice(pre_imputation_df, m = 5, maxit = 50, seed = 826,
             method = "rf", printFlag = TRUE)
}
if (RUN_MICE) {
  message("\n[6/6] 准备运行原始 MICE ...")
  if (!nzchar(PRE_IMPUTATION_FILE) || !file.exists(PRE_IMPUTATION_FILE)) {
    stop("RUN_MICE=TRUE，但 PRE_IMPUTATION_FILE 未设置为有效的插补前数据文件。")
  }
  unimp <- read_analysis_file(PRE_IMPUTATION_FILE)
  imp <- run_original_mice(unimp)
  saveRDS(imp, file.path(OUTPUT_DIR, "mice_object_LOCAL_ONLY.rds"))
  warning("原脚本此前使用 complete(imp,1)；这不是五套插补模型的 Rubin 合并推断。")
} else {
  message("\n[6/6] 已跳过 MICE（本次输入为已插补的 Excel）。")
}

message("\n=================== 全部流程结束 ===================")
message("结果文件夹：", normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE))
print(list.files(OUTPUT_DIR))
message("投稿前需核对插补结果合并、RCS 检验方式和 E-value 计算方法与论文是否一致。")
