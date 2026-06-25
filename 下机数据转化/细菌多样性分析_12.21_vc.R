####准备工作####
#设置读取文件所在的位置（注意window下的斜杠与R相反，需要修改为R背景下的斜杠）
setwd("D:/3_实验数据/")
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