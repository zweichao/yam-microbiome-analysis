############################################################
# SourceTracker 原版 R 包源库分析 + 三张图
#
# 使用的包/代码来源：
#   https://github.com/danknights/sourcetracker
#
# 重要区别：
#   这是 danknights/sourcetracker 原版 R 版 SourceTracker，
#   不是 caporaso-lab/sourcetracker2 命令行版。
#
# 本脚本从三个基础文件开始：
#   feature_table.csv
#   sample_table.csv
#   tax_table.csv
#
# 研究设定：
#   Sources = BS, RS, RP, LP
#   Sinks   = RE, TE, BE, SE, LE
#
# 输出内容：
#   1. 全地区合并版 SourceTracker 结果
#   2. 分地区验证版 SourceTracker 结果
#   3. 三张图，每张图均包含：
#      Part 1. 全地区合并版
#      Part 2. 分地区验证版
#
# 三张图：
#   1. Stacked barplot
#   2. Heatmap
#   3. Sankey / alluvial plot
#
# 当前生态位编号：
#   1 = LP 叶表
#   2 = LE 叶内生
#   3 = SE 茎内生
#   4 = BE 珠芽内生
#   5 = RE 根内生
#   6 = TE 块茎内生
#   7 = RP 根表土
#   8 = RS 根际土
#   9 = BS 非根际土
#
# 展示顺序：
#   Sources: BS, RS, RP, LP, Unknown
#   Sinks:   RE, TE, BE, SE, LE
#
# 注意：
#   SourceTracker 原版要求输入必须是整数 counts，
#   所以本脚本使用 feature_table.csv 原始 reads，
#   不使用相对丰度。
############################################################


# ============================================================
# 0. 清空环境并加载 R 包
# ============================================================

rm(list = ls())

packages <- c(
  "tidyverse",
  "stringr",
  "ggalluvial",
  "scales"
)

new_pkg <- packages[!packages %in% installed.packages()[, "Package"]]

if (length(new_pkg)) {
  install.packages(new_pkg)
}

invisible(lapply(packages, library, character.only = TRUE))


# ============================================================
# 1. 设置路径
# ============================================================

base_dir <- "D:/3_实验数据/3_珠芽微生物/6_不同地区铁棍山药内生菌传播机制/4_微生物数据/2_汇总rawData数据/细菌R - 副本"

work_dir <- file.path(base_dir, "basic_files")

setwd(work_dir)

outdir <- work_dir


# ============================================================
# 2. SourceTracker 原版 R 代码加载
# ============================================================
# 推荐做法：
#   下载 https://github.com/danknights/sourcetracker
#   解压后把整个 sourcetracker 文件夹放在 basic_files 目录下：
#
#   basic_files/sourcetracker/src/SourceTracker.r
#
# 如果你放在别的位置，可以手动修改 sourcetracker_repo。

sourcetracker_repo <- file.path(work_dir, "sourcetracker")

load_sourcetracker_original <- function() {

  if (exists("sourcetracker") && exists("predict.sourcetracker")) {
    message("SourceTracker 函数已存在，跳过加载。")
    return(invisible(TRUE))
  }

  env_path <- Sys.getenv("SOURCETRACKER_PATH")

  candidate_files <- c(
    file.path(sourcetracker_repo, "src", "SourceTracker.r"),
    file.path(work_dir, "src", "SourceTracker.r"),
    file.path(work_dir, "SourceTracker.r"),
    ifelse(env_path != "", file.path(env_path, "src", "SourceTracker.r"), NA_character_),
    ifelse(env_path != "", file.path(env_path, "SourceTracker.r"), NA_character_)
  )

  candidate_files <- candidate_files[!is.na(candidate_files)]

  source_file <- candidate_files[file.exists(candidate_files)][1]

  if (!is.na(source_file)) {
    source(source_file)
    message("已加载 SourceTracker 原版 R 代码：", source_file)
    return(invisible(TRUE))
  }

  # 如果本地没有 SourceTracker.r，则尝试从 GitHub 下载
  message("未在本地找到 SourceTracker.r，尝试从 GitHub 下载。")

  download_to <- file.path(work_dir, "SourceTracker.r")

  try_download <- try(
    download.file(
      url = "https://raw.githubusercontent.com/danknights/sourcetracker/master/src/SourceTracker.r",
      destfile = download_to,
      mode = "wb"
    ),
    silent = TRUE
  )

  if (!inherits(try_download, "try-error") && file.exists(download_to)) {
    source(download_to)
    message("已从 GitHub 下载并加载 SourceTracker.r：", download_to)
    return(invisible(TRUE))
  }

  stop(
    "无法加载 SourceTracker 原版 R 代码。\n",
    "请下载 https://github.com/danknights/sourcetracker，",
    "并把 src/SourceTracker.r 放到 basic_files/sourcetracker/src/SourceTracker.r，",
    "或设置 SOURCETRACKER_PATH。"
  )
}

load_sourcetracker_original()


# ============================================================
# 3. SourceTracker 参数
# ============================================================
# 原版 SourceTracker example.r 中常用：
#   alpha1 <- alpha2 <- 0.001
#
# rarefaction_depth:
#   原函数默认 1000。
#   如果不想稀释，可以改为 NULL。
#   为了与原版默认逻辑接近，这里先用 1000。

alpha1 <- 0.001
alpha2 <- 0.001
beta <- 10

burnin <- 100
nrestarts <- 10
ndraws.per.restart <- 1
delay <- 10

rarefaction_depth <- 1000

set.seed(123)


# ============================================================
# 4. 分组、顺序和颜色
# ============================================================

region_order <- c("BD", "HZ", "WXR", "WXS", "YC", "YJ")

source_order <- c("BS", "RS", "RP", "LP", "Unknown")
source_niches <- c("BS", "RS", "RP", "LP")

sink_order <- c("RE", "TE", "BE", "SE", "LE")

source_colors <- c(
  "BS" = "#888888",
  "RS" = "#946258",
  "RP" = "#9B72C2",
  "LP" = "#5094C3",
  "Unknown" = "#D9D9D9"
)

sink_colors <- c(
  "RE" = "#D83939",
  "TE" = "#FC8F2B",
  "BE" = "#A8E39D",
  "SE" = "#53AA69",
  "LE" = "#48CAD7"
)

source_labels <- c(
  "BS" = "BS\nBulk soil",
  "RS" = "RS\nRhizosphere soil",
  "RP" = "RP\nRoot surface soil",
  "LP" = "LP\nLeaf phyllosphere",
  "Unknown" = "Unknown"
)

sink_labels <- c(
  "RE" = "RE\nRoot endosphere",
  "TE" = "TE\nTuber endosphere",
  "BE" = "BE\nBulbil endosphere",
  "SE" = "SE\nStem endosphere",
  "LE" = "LE\nLeaf endosphere"
)


# ============================================================
# 5. 读取基础数据
# ============================================================

feature_table <- read.csv(
  "feature_table.csv",
  row.names = 1,
  check.names = FALSE
)

sample_table <- read.csv(
  "sample_table.csv",
  row.names = 1,
  check.names = FALSE
)

tax_table <- read.csv(
  "tax_table.csv",
  row.names = 1,
  check.names = FALSE
)

# 确保 feature_table 是整数 counts
feature_table <- round(as.matrix(feature_table))
storage.mode(feature_table) <- "integer"


# ============================================================
# 6. 解析样本信息
# ============================================================

niche_map <- c(
  "1" = "LP",
  "2" = "LE",
  "3" = "SE",
  "4" = "BE",
  "5" = "RE",
  "6" = "TE",
  "7" = "RP",
  "8" = "RS",
  "9" = "BS"
)

parse_region <- function(x) {
  stringr::str_extract(x, "BD|HZ|WXR|WXS|YC|YJ")
}

parse_niche_no <- function(x) {
  parsed <- stringr::str_match(
    x,
    "(BD|HZ|WXR|WXS|YC|YJ)[_\\.-]*([1-9])"
  )
  parsed[, 3]
}

sample_table$SampleID <- rownames(sample_table)

if (!"Group" %in% colnames(sample_table)) {
  sample_table$Group <- sample_table$SampleID
}

sample_table$Region <- parse_region(sample_table$Group)
sample_table$Niche_no <- parse_niche_no(sample_table$Group)
sample_table$Niche <- niche_map[sample_table$Niche_no]

sample_table$Region <- factor(sample_table$Region, levels = region_order)
sample_table$Niche <- factor(
  sample_table$Niche,
  levels = c("LP", "LE", "SE", "BE", "TE", "RE", "RP", "RS", "BS")
)

if (any(is.na(sample_table$Region)) | any(is.na(sample_table$Niche))) {

  failed_sample <- sample_table %>%
    dplyr::filter(is.na(Region) | is.na(Niche)) %>%
    dplyr::select(SampleID, Group, Region, Niche_no, Niche)

  print(failed_sample)
  stop("存在无法解析 Region 或 Niche 的样本，请检查 sample_table.csv 的 Group 列或样本名。")
}


# ============================================================
# 7. 对齐数据并删除零深度样本
# ============================================================

common_samples <- intersect(colnames(feature_table), rownames(sample_table))

if (length(common_samples) == 0) {
  stop("feature_table 的列名与 sample_table 的样本名没有交集，请检查样本名。")
}

feature_table <- feature_table[, common_samples, drop = FALSE]
sample_table <- sample_table[common_samples, , drop = FALSE]

sample_sums <- colSums(feature_table)

if (any(sample_sums == 0)) {
  zero_samples <- names(sample_sums)[sample_sums == 0]
  message("删除总 reads 为 0 的样本：", paste(zero_samples, collapse = ", "))
}

feature_table <- feature_table[, sample_sums > 0, drop = FALSE]
sample_table <- sample_table[colnames(feature_table), , drop = FALSE]

cat("总样本数：", ncol(feature_table), "\n")
cat("总 ASV 数：", nrow(feature_table), "\n")
cat("地区 × 生态位样本数：\n")
print(table(sample_table$Region, sample_table$Niche))


# ============================================================
# 8. SourceTracker 运行函数
# ============================================================
# SourceTracker 原版要求输入是：
#   samples × features 的整数矩阵
# 所以这里把 feature_table 转置。

run_sourcetracker_model <- function(
    model_name,
    sample_ids,
    output_dir
) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  sample_ids <- intersect(sample_ids, colnames(feature_table))

  metadata <- sample_table[sample_ids, , drop = FALSE]

  metadata$SourceSink <- dplyr::case_when(
    as.character(metadata$Niche) %in% source_niches ~ "source",
    as.character(metadata$Niche) %in% sink_order ~ "sink",
    TRUE ~ NA_character_
  )

  metadata$Env <- dplyr::case_when(
    metadata$SourceSink == "source" ~ as.character(metadata$Niche),
    metadata$SourceSink == "sink" ~ as.character(metadata$Niche),
    TRUE ~ NA_character_
  )

  metadata <- metadata %>%
    dplyr::filter(!is.na(SourceSink)) %>%
    dplyr::mutate(
      SourceSink = factor(SourceSink, levels = c("source", "sink")),
      Env = as.character(Env)
    )

  if (nrow(metadata) == 0) {
    warning("模型 ", model_name, " 没有 source/sink 样本，跳过。")
    return(NULL)
  }

  source_ids <- metadata$SampleID[metadata$SourceSink == "source"]
  sink_ids <- metadata$SampleID[metadata$SourceSink == "sink"]

  if (length(source_ids) < 4) {
    warning("模型 ", model_name, " source 样本数过少，跳过。")
    return(NULL)
  }

  if (length(sink_ids) < 1) {
    warning("模型 ", model_name, " 没有 sink 样本，跳过。")
    return(NULL)
  }

  otu_sample_feature <- t(feature_table[, metadata$SampleID, drop = FALSE])

  # 删除所有样本中全 0 的 ASV，减少计算量
  otu_sample_feature <- otu_sample_feature[, colSums(otu_sample_feature) > 0, drop = FALSE]

  train_ix <- which(metadata$SourceSink == "source")
  test_ix <- which(metadata$SourceSink == "sink")

  train_otu <- otu_sample_feature[train_ix, , drop = FALSE]
  test_otu <- otu_sample_feature[test_ix, , drop = FALSE]

  envs_train <- metadata$Env[train_ix]

  cat("\n====================================================\n")
  cat("运行 SourceTracker 模型：", model_name, "\n")
  cat("source 样本数：", nrow(train_otu), "\n")
  cat("sink 样本数：", nrow(test_otu), "\n")
  cat("ASV 数：", ncol(train_otu), "\n")
  cat("source 类型：", paste(sort(unique(envs_train)), collapse = ", "), "\n")
  cat("====================================================\n")

  # 保存输入，方便检查
  write.table(
    otu_sample_feature,
    file = file.path(output_dir, paste0(model_name, "_otu_sample_by_asv_counts.txt")),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )

  write.table(
    metadata,
    file = file.path(output_dir, paste0(model_name, "_metadata.txt")),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )

  # 训练 SourceTracker
  st <- sourcetracker(
    train = train_otu,
    envs = envs_train,
    rarefaction_depth = rarefaction_depth
  )

  # 预测 sink 来源比例
  res <- predict(
    st,
    test = test_otu,
    burnin = burnin,
    nrestarts = nrestarts,
    ndraws.per.restart = ndraws.per.restart,
    delay = delay,
    alpha1 = alpha1,
    alpha2 = alpha2,
    beta = beta,
    rarefaction_depth = rarefaction_depth,
    verbosity = 1,
    full.results = FALSE
  )

  prop <- as.data.frame(res$proportions, check.names = FALSE)
  prop_sd <- as.data.frame(res$proportions_sd, check.names = FALSE)

  # 标准化列名，确保包含 BS, RS, RP, LP, Unknown
  for (s in source_order) {
    if (!s %in% colnames(prop)) {
      prop[[s]] <- 0
    }
    if (!s %in% colnames(prop_sd)) {
      prop_sd[[s]] <- 0
    }
  }

  prop <- prop[, source_order, drop = FALSE]
  prop_sd <- prop_sd[, source_order, drop = FALSE]

  write.table(
    prop,
    file = file.path(output_dir, "mixing_proportions.txt"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )

  write.table(
    prop_sd,
    file = file.path(output_dir, "mixing_proportion_stds.txt"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )

  cat("已输出：", file.path(output_dir, "mixing_proportions.txt"), "\n")

  prop_long <- prop %>%
    tibble::rownames_to_column("SampleID") %>%
    tidyr::pivot_longer(
      cols = all_of(source_order),
      names_to = "Source",
      values_to = "Contribution"
    ) %>%
    dplyr::left_join(
      metadata %>%
        dplyr::select(SampleID, Group, Region, Niche, SourceSink, Env),
      by = "SampleID"
    ) %>%
    dplyr::filter(SourceSink == "sink") %>%
    dplyr::mutate(
      Model = model_name,
      Source = factor(Source, levels = source_order),
      Niche = factor(as.character(Niche), levels = sink_order),
      Region = factor(as.character(Region), levels = region_order)
    )

  write.csv(
    prop_long,
    file = file.path(output_dir, paste0(model_name, "_mixing_proportions_long.csv")),
    row.names = FALSE
  )

  return(prop_long)
}


# ============================================================
# 9. 全地区合并版 SourceTracker
# ============================================================

all_output_dir <- file.path(outdir, "SourceTracker_original_all")

all_sample_ids <- sample_table %>%
  dplyr::filter(Niche %in% c(source_niches, sink_order)) %>%
  dplyr::pull(SampleID)

mix_all_long <- run_sourcetracker_model(
  model_name = "all_region",
  sample_ids = all_sample_ids,
  output_dir = all_output_dir
)


# ============================================================
# 10. 分地区验证版 SourceTracker
# ============================================================

by_region_output_dir <- file.path(outdir, "SourceTracker_original_by_region")

mix_region_long <- purrr::map_dfr(
  region_order,
  function(region_i) {

    sample_ids_i <- sample_table %>%
      dplyr::filter(
        Region == region_i,
        Niche %in% c(source_niches, sink_order)
      ) %>%
      dplyr::pull(SampleID)

    run_sourcetracker_model(
      model_name = paste0("region_", region_i),
      sample_ids = sample_ids_i,
      output_dir = file.path(by_region_output_dir, region_i)
    ) %>%
      dplyr::mutate(RegionModel = region_i)
  }
)

if (!is.null(mix_all_long)) {
  write.csv(
    mix_all_long,
    file.path(outdir, "SourceTracker_original_all_region_long_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.csv"),
    row.names = FALSE
  )
}

if (!is.null(mix_region_long) && nrow(mix_region_long) > 0) {
  mix_region_long$RegionModel <- factor(mix_region_long$RegionModel, levels = region_order)

  write.csv(
    mix_region_long,
    file.path(outdir, "SourceTracker_original_by_region_long_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.csv"),
    row.names = FALSE
  )
}


# ============================================================
# 11. 汇总函数
# ============================================================

make_summary <- function(df, group_vars) {

  df %>%
    dplyr::group_by(dplyr::across(all_of(c(group_vars, "Source")))) %>%
    dplyr::summarise(
      MeanContribution = mean(Contribution, na.rm = TRUE),
      SD = sd(Contribution, na.rm = TRUE),
      SE = SD / sqrt(dplyr::n()),
      N = dplyr::n(),
      .groups = "drop"
    )
}

summary_all <- make_summary(
  mix_all_long,
  group_vars = "Niche"
) %>%
  tidyr::complete(
    Niche = factor(sink_order, levels = sink_order),
    Source = factor(source_order, levels = source_order),
    fill = list(
      MeanContribution = 0,
      SD = 0,
      SE = 0,
      N = 0
    )
  ) %>%
  dplyr::group_by(Niche) %>%
  dplyr::mutate(
    Contribution_plot = MeanContribution / sum(MeanContribution, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Contribution_plot = ifelse(is.na(Contribution_plot), 0, Contribution_plot),
    Niche = factor(as.character(Niche), levels = sink_order),
    Source = factor(as.character(Source), levels = source_order)
  )

summary_region <- make_summary(
  mix_region_long,
  group_vars = c("RegionModel", "Niche")
) %>%
  tidyr::complete(
    RegionModel = factor(region_order, levels = region_order),
    Niche = factor(sink_order, levels = sink_order),
    Source = factor(source_order, levels = source_order),
    fill = list(
      MeanContribution = 0,
      SD = 0,
      SE = 0,
      N = 0
    )
  ) %>%
  dplyr::group_by(RegionModel, Niche) %>%
  dplyr::mutate(
    Contribution_plot = MeanContribution / sum(MeanContribution, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Contribution_plot = ifelse(is.na(Contribution_plot), 0, Contribution_plot),
    RegionModel = factor(as.character(RegionModel), levels = region_order),
    Niche = factor(as.character(Niche), levels = sink_order),
    Source = factor(as.character(Source), levels = source_order)
  )

write.csv(
  summary_all,
  file.path(outdir, "SourceTracker_original_all_region_summary_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.csv"),
  row.names = FALSE
)

write.csv(
  summary_region,
  file.path(outdir, "SourceTracker_original_by_region_summary_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.csv"),
  row.names = FALSE
)


# ============================================================
# 12. 图 1：Stacked barplot
# ============================================================

p_bar_all <- ggplot(
  summary_all,
  aes(x = Niche, y = Contribution_plot, fill = Source)
) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = source_colors, breaks = source_order, labels = source_labels) +
  scale_x_discrete(labels = sink_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Estimated source contribution", fill = "Source") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

# 预览图片
print(p_bar_all)

ggsave(file.path(outdir, "SourceTracker_original_stacked_barplot_all_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.pdf"), p_bar_all, width = 7.5, height = 4.8)
ggsave(file.path(outdir, "SourceTracker_original_stacked_barplot_all_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.png"), p_bar_all, width = 7.5, height = 4.8, dpi = 600)
ggsave(file.path(outdir, "SourceTracker_original_stacked_barplot_all_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.svg"), p_bar_all, width = 7.5, height = 4.8)

p_bar_region <- ggplot(
  summary_region,
  aes(x = Niche, y = Contribution_plot, fill = Source)
) +
  geom_col(width = 0.72, color = "white", linewidth = 0.20) +
  facet_wrap(~ RegionModel, nrow = 2) +
  scale_fill_manual(values = source_colors, breaks = source_order, labels = source_labels) +
  scale_x_discrete(labels = sink_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Estimated source contribution", fill = "Source") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.y = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.position = "right"
  )

# 预览图片
print(p_bar_region)

ggsave(file.path(outdir, "SourceTracker_original_stacked_barplot_by_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.pdf"), p_bar_region, width = 10.5, height = 7.2)
ggsave(file.path(outdir, "SourceTracker_original_stacked_barplot_by_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.png"), p_bar_region, width = 10.5, height = 7.2, dpi = 600)
ggsave(file.path(outdir, "SourceTracker_original_stacked_barplot_by_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.svg"), p_bar_region, width = 10.5, height = 7.2)


# ============================================================
# 13. 图 2：Heatmap
# ============================================================

p_heat_all <- ggplot(
  summary_all,
  aes(x = Source, y = Niche, fill = Contribution_plot)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = scales::percent(Contribution_plot, accuracy = 0.1)), size = 3.4, color = "black") +
  scale_x_discrete(labels = source_labels, position = "top") +
  scale_y_discrete(labels = sink_labels) +
  scale_fill_gradient(low = "white", high = "#D83939", labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(x = NULL, y = NULL, fill = "Contribution") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

# 预览图片
print(p_heat_all)

ggsave(file.path(outdir, "SourceTracker_original_heatmap_all_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.pdf"), p_heat_all, width = 7.2, height = 4.8)
ggsave(file.path(outdir, "SourceTracker_original_heatmap_all_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.png"), p_heat_all, width = 7.2, height = 4.8, dpi = 600)
ggsave(file.path(outdir, "SourceTracker_original_heatmap_all_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.svg"), p_heat_all, width = 7.2, height = 4.8)

p_heat_region <- ggplot(
  summary_region,
  aes(x = Source, y = Niche, fill = Contribution_plot)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = scales::percent(Contribution_plot, accuracy = 0.1)), size = 2.4, color = "black") +
  facet_wrap(~ RegionModel, nrow = 2) +
  scale_x_discrete(
    labels = c("BS" = "BS", "RS" = "RS", "RP" = "RP", "LP" = "LP", "Unknown" = "Unknown"),
    position = "top"
  ) +
  scale_y_discrete(labels = c("RE" = "RE", "TE" = "TE", "BE" = "BE", "SE" = "SE", "LE" = "LE")) +
  scale_fill_gradient(low = "white", high = "#D83939", labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(x = NULL, y = NULL, fill = "Contribution") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 0),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.position = "right"
  )

# 预览图片
print(p_heat_region)

ggsave(file.path(outdir, "SourceTracker_original_heatmap_by_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.pdf"), p_heat_region, width = 10.5, height = 7.2)
ggsave(file.path(outdir, "SourceTracker_original_heatmap_by_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.png"), p_heat_region, width = 10.5, height = 7.2, dpi = 600)
ggsave(file.path(outdir, "SourceTracker_original_heatmap_by_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.svg"), p_heat_region, width = 10.5, height = 7.2)


# ============================================================
# 14. 图 3：Sankey / alluvial
# ============================================================

sankey_all <- summary_all %>%
  dplyr::filter(Contribution_plot > 0) %>%
  dplyr::mutate(
    Source = factor(as.character(Source), levels = source_order),
    Niche = factor(as.character(Niche), levels = sink_order)
  )

p_sankey_all <- ggplot(
  sankey_all,
  aes(axis1 = Source, axis2 = Niche, y = Contribution_plot)
) +
  ggalluvial::geom_alluvium(
    aes(fill = Source),
    width = 1 / 8,
    alpha = 0.78,
    color = "white",
    linewidth = 0.15
  ) +
  ggalluvial::geom_stratum(
    width = 1 / 8,
    fill = "grey95",
    color = "grey40",
    linewidth = 0.35
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 3.5,
    color = "black"
  ) +
  scale_x_discrete(limits = c("Source", "Sink"), expand = c(0.12, 0.08)) +
  scale_fill_manual(values = source_colors, breaks = source_order, labels = source_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, y = "Estimated source contribution", fill = "Source") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

# 预览图片
print(p_sankey_all)

ggsave(file.path(outdir, "SourceTracker_original_sankey_all_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.pdf"), p_sankey_all, width = 7.8, height = 5.2)
ggsave(file.path(outdir, "SourceTracker_original_sankey_all_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.png"), p_sankey_all, width = 7.8, height = 5.2, dpi = 600)
ggsave(file.path(outdir, "SourceTracker_original_sankey_all_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.svg"), p_sankey_all, width = 7.8, height = 5.2)

sankey_region <- summary_region %>%
  dplyr::filter(Contribution_plot > 0) %>%
  dplyr::mutate(
    RegionModel = factor(as.character(RegionModel), levels = region_order),
    Source = factor(as.character(Source), levels = source_order),
    Niche = factor(as.character(Niche), levels = sink_order)
  )

p_sankey_region <- ggplot(
  sankey_region,
  aes(axis1 = Source, axis2 = Niche, y = Contribution_plot)
) +
  ggalluvial::geom_alluvium(
    aes(fill = Source),
    width = 1 / 8,
    alpha = 0.78,
    color = "white",
    linewidth = 0.10
  ) +
  ggalluvial::geom_stratum(
    width = 1 / 8,
    fill = "grey95",
    color = "grey40",
    linewidth = 0.25
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 2.7,
    color = "black"
  ) +
  facet_wrap(~ RegionModel, nrow = 2) +
  scale_x_discrete(limits = c("Source", "Sink"), expand = c(0.12, 0.08)) +
  scale_fill_manual(values = source_colors, breaks = source_order, labels = source_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, y = "Estimated source contribution", fill = "Source") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, color = "black", face = "bold"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.position = "right"
  )

# 预览图片
print(p_sankey_region)

ggsave(file.path(outdir, "SourceTracker_original_sankey_by_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.pdf"), p_sankey_region, width = 11, height = 7.5)
ggsave(file.path(outdir, "SourceTracker_original_sankey_by_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.png"), p_sankey_region, width = 11, height = 7.5, dpi = 600)
ggsave(file.path(outdir, "SourceTracker_original_sankey_by_region_BS_RS_RP_LP_to_RE_TE_BE_SE_LE.svg"), p_sankey_region, width = 11, height = 7.5)


# ============================================================
# 15. 完成提示
# ============================================================

message("\nSourceTracker 原版 R 源库分析与三张图全部完成。")
message("输出目录：", outdir)
message("全地区结果目录：", all_output_dir)
message("分地区结果目录：", by_region_output_dir)
message("Source = BS, RS, RP, LP; Sink = RE, TE, BE, SE, LE")
