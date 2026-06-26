############################################################
# 铁棍山药微生物群落分析
# 9 个生态位 ASV 共现网络 — ggClusterNet-style Gephi 输出版
#
# 参照方法描述：
#   Correlation networks were constructed using ggClusterNet-style workflow.
#   Spearman correlation was chosen because it does not rely on normality
#   or linearity assumptions.
#   Correlations with |Spearman r| > 0.6 and P < 0.05 were retained.
#   Network visualization was performed with Gephi.
#
# 本脚本说明：
#   1. 以 ASV 为节点；
#   2. 在每个生态位内部独立计算 ASV-ASV Spearman 相关；
#   3. 保留 |r| > 0.6 且 P < 0.05 的边；
#   4. 删除孤立节点；
#   5. 只导出 Gephi 可打开的 graphml 文件；
#   6. 不生成 R 网络图、节点表和边表。
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
# 输出顺序：
#   LP, LE, SE, BE, TE, RE, RP, RS, BS
#
# 主要输出：
#   ggClusterNet_style_network_LP_spearman_r06_p005_updated_5RE_6TE.graphml
#   ...
#   ggClusterNet_style_network_BS_spearman_r06_p005_updated_5RE_6TE.graphml
############################################################


# ============================================================
# 0. 清空环境并加载 R 包
# ============================================================

rm(list = ls())

packages <- c(
  "tidyverse",
  "stringr",
  "psych",
  "igraph",
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
# 2. 网络构建参数
# ============================================================

# 与参考方法一致：Spearman |r| > 0.6，P < 0.05
cor_method <- "spearman"
r_threshold <- 0.6
p_threshold <- 0.05

# 参考方法中写的是 P < 0.05，未说明 FDR；
# 因此这里默认不做 FDR 校正，使用原始 P 值。
# 如果后续需要更严格版本，可改为 "fdr"。
p_adjust_method <- "none"

# 低丰度 ASV 过滤：用于减少稀有 ASV 噪音。
# 如果希望完全按参考文字不加丰度过滤，可改为 0。
min_mean_rel_abund <- 0.0005

# 出现率过滤：当前默认不启用。
# 如果想至少出现在 20% 样本中，可改为 0.20。
prevalence_cutoff <- 0

# presence 判定阈值
detection_cutoff <- 0

# 删除没有显著边的孤立节点
delete_single <- TRUE

# 是否标准化 ASV 丰度
# Spearman 基于秩相关，标准化对相关方向影响通常不大；
# 默认 FALSE，更直接对应相对丰度数据。
scale_data <- FALSE

set.seed(123)


# ============================================================
# 3. 地区和生态位信息
# ============================================================

region_order <- c("BD", "HZ", "WXR", "WXS", "YC", "YJ")

# 关键：5 = RE，6 = TE
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

niche_order <- c("LP", "LE", "SE", "BE", "TE", "RE", "RP", "RS", "BS")

niche_full <- c(
  "LP" = "Leaf phyllosphere",
  "LE" = "Leaf endosphere",
  "SE" = "Stem endosphere",
  "BE" = "Bulbil endosphere",
  "TE" = "Tuber endosphere",
  "RE" = "Root endosphere",
  "RP" = "Root surface soil",
  "RS" = "Rhizosphere soil",
  "BS" = "Bulk soil"
)


# ============================================================
# 4. 读取数据
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


# ============================================================
# 5. 样本信息解析
# ============================================================

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
sample_table$Niche <- factor(sample_table$Niche, levels = niche_order)

if (any(is.na(sample_table$Region)) | any(is.na(sample_table$Niche))) {
  
  failed_sample <- sample_table %>%
    dplyr::filter(is.na(Region) | is.na(Niche)) %>%
    dplyr::select(SampleID, Group, Region, Niche_no, Niche)
  
  print(failed_sample)
  stop("存在无法解析 Region 或 Niche 的样本，请检查 sample_table.csv 的 Group 列或样本名。")
}


# ============================================================
# 6. 对齐 feature_table、sample_table 和 tax_table
# ============================================================

common_samples <- intersect(colnames(feature_table), rownames(sample_table))

if (length(common_samples) == 0) {
  stop("feature_table 的列名与 sample_table 的样本名没有交集，请检查样本名。")
}

feature_table <- feature_table[, common_samples, drop = FALSE]
sample_table <- sample_table[common_samples, , drop = FALSE]

common_asv <- intersect(rownames(feature_table), rownames(tax_table))

if (length(common_asv) == 0) {
  stop("feature_table 的行名与 tax_table 的行名没有交集，请检查 ASV/OTU 名称。")
}

feature_table <- feature_table[common_asv, , drop = FALSE]
tax_table <- tax_table[common_asv, , drop = FALSE]

# 删除总 reads 为 0 的样本
sample_sums <- colSums(feature_table)

if (any(sample_sums == 0)) {
  zero_samples <- names(sample_sums)[sample_sums == 0]
  message("删除总丰度为 0 的样本：", paste(zero_samples, collapse = ", "))
}

feature_table <- feature_table[, sample_sums > 0, drop = FALSE]
sample_table <- sample_table[colnames(feature_table), , drop = FALSE]
feature_table <- feature_table[, rownames(sample_table), drop = FALSE]

# 样本内相对丰度
feature_rel <- sweep(feature_table, 2, colSums(feature_table), "/")
feature_rel[is.na(feature_rel)] <- 0

cat("总样本数：", ncol(feature_table), "\n")
cat("总 ASV 数：", nrow(feature_table), "\n")
cat("地区 × 生态位样本数：\n")
print(table(sample_table$Region, sample_table$Niche))


# ============================================================
# 7. 分类注释整理
# ============================================================

clean_tax_value <- function(x) {
  x <- as.character(x)
  x <- stringr::str_remove(x, "^[a-z]__")
  x[x == ""] <- NA
  x
}

get_tax_col <- function(tax_df, candidates) {
  x <- intersect(candidates, colnames(tax_df))
  if (length(x) == 0) {
    return(rep(NA_character_, nrow(tax_df)))
  } else {
    return(tax_df[[x[1]]])
  }
}

make_tax_info <- function(tax_df) {
  
  tax_info <- tax_df %>%
    tibble::rownames_to_column("ASV")
  
  tax_info$Kingdom <- clean_tax_value(get_tax_col(tax_df, c("Kingdom", "kingdom", "Domain", "domain")))
  tax_info$Phylum  <- clean_tax_value(get_tax_col(tax_df, c("Phylum", "phylum")))
  tax_info$Class   <- clean_tax_value(get_tax_col(tax_df, c("Class", "class")))
  tax_info$Order   <- clean_tax_value(get_tax_col(tax_df, c("Order", "order")))
  tax_info$Family  <- clean_tax_value(get_tax_col(tax_df, c("Family", "family", "familiy")))
  tax_info$Genus   <- clean_tax_value(get_tax_col(tax_df, c("Genus", "genus")))
  tax_info$Species <- clean_tax_value(get_tax_col(tax_df, c("Species", "species", "specise")))
  
  tax_info <- tax_info %>%
    mutate(
      Taxon_label = case_when(
        !is.na(Genus) & Genus != "" ~ paste0("g__", Genus),
        !is.na(Family) & Family != "" ~ paste0("f__", Family),
        !is.na(Order) & Order != "" ~ paste0("o__", Order),
        !is.na(Class) & Class != "" ~ paste0("c__", Class),
        !is.na(Phylum) & Phylum != "" ~ paste0("p__", Phylum),
        TRUE ~ "Unclassified"
      )
    )
  
  return(tax_info)
}

tax_info <- make_tax_info(tax_table)


# ============================================================
# 8. GraphML 属性清理函数
# ============================================================

sanitize_graph_for_graphml <- function(g) {
  
  v_attrs <- igraph::vertex_attr_names(g)
  
  for (attr_i in v_attrs) {
    
    value_i <- igraph::vertex_attr(g, attr_i)
    
    if (is.logical(value_i)) {
      value_i <- as.character(value_i)
      value_i[is.na(value_i)] <- "FALSE"
    } else if (is.numeric(value_i) | is.integer(value_i)) {
      value_i[is.na(value_i)] <- 0
    } else {
      value_i <- as.character(value_i)
      value_i[is.na(value_i)] <- ""
    }
    
    g <- igraph::set_vertex_attr(g, attr_i, value = value_i)
  }
  
  e_attrs <- igraph::edge_attr_names(g)
  
  for (attr_i in e_attrs) {
    
    value_i <- igraph::edge_attr(g, attr_i)
    
    if (is.logical(value_i)) {
      value_i <- as.character(value_i)
      value_i[is.na(value_i)] <- "FALSE"
    } else if (is.numeric(value_i) | is.integer(value_i)) {
      value_i[is.na(value_i)] <- 0
    } else {
      value_i <- as.character(value_i)
      value_i[is.na(value_i)] <- ""
    }
    
    g <- igraph::set_edge_attr(g, attr_i, value = value_i)
  }
  
  return(g)
}


# ============================================================
# 9. 单个生态位网络构建函数
# ============================================================

build_one_niche_network <- function(niche_name) {
  
  cat("\n====================================================\n")
  cat("Processing niche:", niche_name, "\n")
  cat("====================================================\n")
  
  sample_ids <- sample_table %>%
    dplyr::filter(Niche == niche_name) %>%
    dplyr::pull(SampleID)
  
  sample_ids <- intersect(sample_ids, colnames(feature_rel))
  
  if (length(sample_ids) < 4) {
    warning(paste0("生态位 ", niche_name, " 样本数少于 4，跳过。"))
    return(NULL)
  }
  
  # ASV × sample 相对丰度表
  rel_sub <- feature_rel[, sample_ids, drop = FALSE]
  
  # 过滤全 0 ASV
  rel_sub <- rel_sub[rowSums(rel_sub) > 0, , drop = FALSE]
  
  # 平均相对丰度过滤
  mean_rel <- rowMeans(rel_sub, na.rm = TRUE)
  prevalence <- rowMeans(rel_sub > detection_cutoff, na.rm = TRUE)
  
  keep_asv <- names(mean_rel)[
    mean_rel > min_mean_rel_abund &
      prevalence >= prevalence_cutoff
  ]
  
  if (length(keep_asv) < 3) {
    warning(paste0("生态位 ", niche_name, " 丰度/出现率过滤后 ASV 少于 3 个，跳过。"))
    return(NULL)
  }
  
  rel_sub <- rel_sub[keep_asv, , drop = FALSE]
  mean_rel <- mean_rel[keep_asv]
  prevalence <- prevalence[keep_asv]
  
  # 转为 样本 × ASV
  spe <- t(rel_sub) %>%
    as.data.frame(check.names = FALSE)
  
  # 去除标准差为 0 或有效样本数不足的 ASV
  good_asv <- apply(
    spe,
    2,
    function(x) {
      sd(x, na.rm = TRUE) > 1e-8 && sum(!is.na(x)) >= 3
    }
  )
  
  if (sum(!good_asv) > 0) {
    cat("去除零方差或有效样本不足 ASV 数：", sum(!good_asv), "\n")
  }
  
  spe <- spe[, good_asv, drop = FALSE]
  
  if (ncol(spe) < 3) {
    warning(paste0("生态位 ", niche_name, " 去除零方差后 ASV 少于 3 个，跳过。"))
    return(NULL)
  }
  
  if (scale_data) {
    spe <- scale(spe) %>%
      as.data.frame(check.names = FALSE)
  }
  
  cat("用于网络分析的样本数：", nrow(spe), "\n")
  cat("用于网络分析的 ASV 数：", ncol(spe), "\n")
  
  # 计算 Spearman 相关
  # adjust = "none" 保持和参考方法中 P < 0.05 一致
  cor_res <- psych::corr.test(
    spe,
    method = cor_method,
    adjust = "none"
  )
  
  r_mat <- cor_res$r
  p_mat <- cor_res$p
  
  # 可选 p 值校正
  p_use <- p_mat
  
  upper_idx <- upper.tri(p_mat)
  
  if (!is.null(p_adjust_method) && p_adjust_method != "none") {
    
    p_adj_vec <- p.adjust(
      p_mat[upper_idx],
      method = p_adjust_method
    )
    
    p_use[upper_idx] <- p_adj_vec
    p_use[lower.tri(p_use)] <- t(p_use)[lower.tri(p_use)]
  }
  
  # 按 |r| > 0.6 和 P < 0.05 筛选
  r_filter <- r_mat
  
  r_filter[
    abs(r_filter) <= r_threshold |
      p_use >= p_threshold
  ] <- 0
  
  diag(r_filter) <- 0
  
  # 构建网络
  g <- igraph::graph_from_adjacency_matrix(
    r_filter,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  if (delete_single) {
    g <- igraph::delete_vertices(
      g,
      names(igraph::degree(g)[igraph::degree(g) == 0])
    )
  }
  
  if (igraph::vcount(g) == 0 | igraph::ecount(g) == 0) {
    warning(paste0("生态位 ", niche_name, " 在当前阈值下没有显著边，未导出 graphml。"))
    return(NULL)
  }
  
  # 边属性：保留原始正负相关；weight 用绝对值，方便 Gephi 调边宽
  igraph::E(g)$Correlation <- igraph::E(g)$weight
  igraph::E(g)$R <- ifelse(igraph::E(g)$Correlation > 0, "pos", "neg")
  igraph::E(g)$Interaction <- ifelse(igraph::E(g)$Correlation > 0, "Positive", "Negative")
  igraph::E(g)$Edge_color <- ifelse(igraph::E(g)$Correlation > 0, "#D95F8D", "#4D6CFA")
  igraph::E(g)$Weight <- abs(igraph::E(g)$Correlation)
  igraph::E(g)$weight <- abs(igraph::E(g)$Correlation)
  
  # 节点属性
  node_tax <- tax_info %>%
    dplyr::filter(ASV %in% igraph::V(g)$name) %>%
    dplyr::arrange(match(ASV, igraph::V(g)$name))
  
  igraph::V(g)$ASV <- igraph::V(g)$name
  igraph::V(g)$Niche <- niche_name
  igraph::V(g)$Niche_full <- niche_full[niche_name]
  
  igraph::V(g)$Kingdom <- node_tax$Kingdom
  igraph::V(g)$Phylum <- node_tax$Phylum
  igraph::V(g)$Class <- node_tax$Class
  igraph::V(g)$Order <- node_tax$Order
  igraph::V(g)$Family <- node_tax$Family
  igraph::V(g)$Genus <- node_tax$Genus
  igraph::V(g)$Species <- node_tax$Species
  igraph::V(g)$Taxon_label <- node_tax$Taxon_label
  
  mean_map <- mean_rel
  prev_map <- prevalence
  
  igraph::V(g)$MeanRelAbund <- mean_map[igraph::V(g)$name]
  igraph::V(g)$Prevalence <- prev_map[igraph::V(g)$name]
  
  igraph::V(g)$Degree <- igraph::degree(g)
  igraph::V(g)$Betweenness <- igraph::betweenness(g, directed = FALSE)
  igraph::V(g)$Closeness <- igraph::closeness(g, normalized = TRUE)
  
  igraph::V(g)$Label <- ifelse(
    !is.na(igraph::V(g)$Taxon_label) & igraph::V(g)$Taxon_label != "",
    paste0(igraph::V(g)$Taxon_label, " | ", igraph::V(g)$ASV),
    igraph::V(g)$ASV
  )
  
  # 小写 label 更容易被 Gephi 识别为标签
  igraph::V(g)$label <- ifelse(
    !is.na(igraph::V(g)$Taxon_label) & igraph::V(g)$Taxon_label != "",
    igraph::V(g)$Taxon_label,
    igraph::V(g)$ASV
  )
  
  igraph::V(g)$Node_size <- ifelse(
    length(unique(igraph::V(g)$Degree)) <= 1,
    30,
    scales::rescale(igraph::V(g)$Degree, to = c(10, 60))
  )
  
  igraph::V(g)$Node_color_by_degree <- igraph::V(g)$Degree
  
  g <- sanitize_graph_for_graphml(g)
  
  # 导出 graphml
  out_file <- file.path(
    outdir,
    paste0(
      "ggClusterNet_style_network_",
      niche_name,
      "_spearman_r06_p005_updated_5RE_6TE.graphml"
    )
  )
  
  igraph::write_graph(
    graph = g,
    file = out_file,
    format = "graphml"
  )
  
  cat("正相关边数：", sum(igraph::E(g)$R == "pos"), "\n")
  cat("负相关边数：", sum(igraph::E(g)$R == "neg"), "\n")
  cat("已导出：", out_file, "\n")
  
  return(g)
}


# ============================================================
# 10. 批量运行 9 个生态位
# ============================================================

network_list <- list()

for (niche_i in niche_order) {
  network_list[[niche_i]] <- build_one_niche_network(niche_i)
}


# ============================================================
# 11. 完成提示
# ============================================================

message("\n9 个生态位 ggClusterNet-style Spearman 网络 graphml 文件导出完成。")
message("输出目录：", outdir)
message("输出文件命名格式：")
message("ggClusterNet_style_network_生态位_spearman_r06_p005_updated_5RE_6TE.graphml")
message("建网参数：")
message("cor_method = ", cor_method)
message("|r| > ", r_threshold)
message("P < ", p_threshold)
message("p_adjust_method = ", p_adjust_method)
message("min_mean_rel_abund = ", min_mean_rel_abund)
message("prevalence_cutoff = ", prevalence_cutoff)
message("delete_single = ", delete_single)
