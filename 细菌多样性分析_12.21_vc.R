####准备工作####
#设置读取文件所在的位置（注意window下的斜杠与R相反，需要修改为R背景下的斜杠）
setwd("D:/3_实验数据/其他数据/1_林梦丽/微生物/细菌/R")
#调用后续所需要的R包，没有的话安装即可，注意不同的包安装的方式有所不一样
library(file2meco)
library(magrittr)
#创建文件夹用于后续的存图及文件
if (!dir.create("Files")) {
  dir.create("Files")
}

if (!dir.create("sample_net_file")) {
  dir.create("sample_net_file")
}
if (!dir.create("Images")) {
  dir.create("Images")
}
####数据前处理####
#导入数据
dataset<-qiime2meco(
  'table.qza',
  sample_table  = "metadata.txt",
  taxonomy_table = 'taxonomy.qza',
  phylo_tree = "rooted-tree.qza",
  rep_fasta = "rep-seqs.qza",auto_tidy=T)
#set.seed 用于固定随机数生成机制，确保实验结果具有可重复性,服务于后续的抽平步骤
set.seed(123)
#设置出图类型与风格
library(ggplot2)
theme_set(theme_classic())
#数据清洗（只保留细菌和古菌）
dataset$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
dataset
#数据清洗（去除叶绿体和线粒体的数据）
dataset$filter_pollution(taxa = c("mitochondria", "chloroplast"))
dataset
dataset$rep_fasta@ranges@width%>%summary()
#可选步骤1：筛选长度大于350bp的序列
dataset$otu_table<-dataset$otu_table[dataset$rep_fasta@ranges@width>350,]
#可选步骤2：筛选在五个以上样本中出现的asv
dataset$otu_table<-dataset$otu_table[dataset$otu_table%>%rowSums(dataset$otu_table>0)>5,]
#完成数据的清洗和筛选后对数据进行整理
dataset$tidy_dataset()
dataset
#抽平
a<-dataset$sample_sums() %>% range
dataset$rarefy_samples(sample.size = a[1])
dataset$sample_sums() %>% range
dataset
rm(a)
#将结果保存
dataset$cal_abund()
dataset$save_table(dirpath = "basic_files", sep = ",")
#将各级别的丰度表写出
dataset$save_abund(dirpath = "taxa_abund")
dataset$save_abund(merge_all = TRUE, sep = "\t", quote = FALSE)
dataset$save_abund(merge_all = TRUE, sep = "\t", rm_un = TRUE, rm_pattern = "__$|Sedis$", quote = FALSE)
#alpha和beta多样性计算
dataset$cal_alphadiv(PD = TRUE)
dataset$save_alphadiv(dirpath = "alpha_diversity")#将结果存储
invisible(dataset$cal_betadiv(unifrac = FALSE))
dataset$cal_betadiv(unifrac = TRUE)
dataset$save_betadiv(dirpath = "beta_diversity")
####物种组成作图####
#前10个门水平堆叠图
library(microeco)
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8)
g1<-t1$plot_bar(others_color = "grey70", facet = "Group", xtext_keep = TRUE, 
                legend_text_italic = FALSE,order_x = row.names(dataset$sample_table),
                xtext_size =7)
g1 +  theme(axis.text.x = element_text(angle = 45))
ggsave("Images/plot_bar.pdf",width = 8,height = 5)
#按组前10个门平均
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
g1 + theme_classic() + theme(axis.title.y = element_text(size = 18))
ggsave("Images/plot_Phylum_bar_mean.pdf",width = 5,height = 8)
#前20个属堆叠图
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 20)
g1<-t1$plot_bar(color_values = c(RColorBrewer::brewer.pal(12, "Paired"),
                                 RColorBrewer::brewer.pal(12, "Set3")),
                others_color = "grey70", facet = "Group", xtext_keep = TRUE, 
                legend_text_italic = FALSE,xtext_size =7)
g1 +  theme(axis.text.x = element_text(angle = 45))+guides(fill=guide_legend(reverse=T,title=NULL,ncol=1))
ggsave("Images/plot_Genus_bar.pdf",width = 12,height = 6)
#按组平均的属水平堆叠图
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 20, groupmean = "Group")
g1 <- t1$plot_bar(color_values = c(RColorBrewer::brewer.pal(12, "Paired"),RColorBrewer::brewer.pal(12, "Set3")),others_color = "grey70", legend_text_italic = FALSE)
g1 + theme_classic() + theme(axis.title.y = element_text(size = 18))
ggsave("Images/plot_Genus_bar_mean.pdf",width = 8,height = 5)
#前15个属在各组的相对丰度分布
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 15)
t1$plot_box(group = "Group")
ggsave("Images/plot_Genus_box.pdf",width = 8,height = 5)
#前40个属的热图
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 40)
t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)
ggsave("Images/plot_heatmap.pdf",width = 10,height = 6)
#前8个属在样品中的分布情况
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 8)
t1$plot_line()
ggsave("Images/plot_sampleline_type.pdf",width = 19.2,height = 10.8)
#前8个属在各组的分布情况
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 8, group = "Group")
t1$plot_line(position = position_dodge(0.3))
ggsave("Images/plot_line_Group.pdf",width = 8,height = 5)
#前6个门在各组分布的饼图
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 6, groupmean = "Group")
t1$plot_pie(facet_nrow = 1)
ggsave("Images/plot_pie.pdf",width = 8,height = 5)
#韦恩图
dataset1 <- dataset$merge_samples(group = "Group")
t1 <- trans_venn$new(dataset1, ratio = NULL)
if (t1[["colnumber"]]<=5){
  t1$plot_venn()
  ggsave("Images/plot_venn5.pdf",width = 8,height = 5)
  
  }else{
  t1 <- trans_venn$new(dataset1)
  t1$plot_venn(petal_plot = TRUE, petal_center_size = 50, petal_r =1.5, 
               petal_a = 3, petal_move_xy = 3.8, petal_color_center = "#BEBADA")
  ggsave("Images/trans_venn_2.pdf",width = 8,height = 5)
  }
###组间差异分析####
#alpha多样性
t1$cal_diff(method = "KW")#整体差异
head(t1$res_diff)
suppressWarnings(t1$cal_diff(method = "KW"))
pander::pander(head(t1$res_diff[, -3]))

t1$cal_diff(method = "KW_dunn", KW_dunn_letter = FALSE)
head(t1$res_diff)
t1$res_diff%<>% base::subset(Significance != "ns")

a<-unique(t1$res_diff$Measure)%>%.[-which(df[3]==1)]#秩和检验作图
for (i in a) {
  t1$plot_alpha(measure = i, xtext_size = 15)
  filename=paste0("Images/plot_alpha_kw_",i,".pdf")
  ggsave(file=filename,width = 5,height = 8)
  
}
#beta多样性
t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = "bray")#PCoA,PCA和NMDS
a<-c("PCoA","PCA","NMDS")
for (i in a) {
  t1$cal_ordination(ordination = i)
  # t1$res_ordination is the ordination result list
  t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))
  filename=paste0("Images/plot_ordination_",i,".pdf")
  ggsave(file=filename,width = 5,height = 8)
  
}
####LEfSe差异分析####
t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Group", 
                     alpha = 0.05, lefse_subgroup = NULL)
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
t1$plot_diff_bar(threshold = 2,color_values=c("#4DBBD5FF", "#00A087FF", "#3C5488FF", "#8491B4FF", "#91D1C2FF" ))
ggsave("Images/plot_lefse_threshold2_bar1.pdf",width = 19.2,height = 10.8)
t1$res_diff %<>% base::subset(LDA>2)
write.csv(t1$res_diff,"Files/Lefse_taxa_all.csv")

# we show 30 taxa with the highest LDA (log10)
t1$plot_diff_bar(use_number = 1:50, width = 0.8)
#ggsave("Images/plot_lefse_log10_bar.pdf",width = 8,height = 5)
# show part of the table
#t1$res_diff[1:5, c(1, 3, 4, 6)]
#write.csv(t1$res_diff,"Files/res_lefse.csv")
#res_lefse <- read.csv("Files/res_lefse.csv")
#pander::pander(res_lefse[1:5, c(1, 3, 4, 6)])
####7.2 绘制LEfSe检测到的标志物的丰度####
#Then, we plot the abundance of biomarkers detected by LEfSe.
t1$plot_diff_abund(use_number = 1:50)
ggsave("Images/plot_lefse_diff_abund.pdf",width = 10,height = 5)


