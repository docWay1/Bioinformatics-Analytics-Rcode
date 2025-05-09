#2024-04-27开始重新学习单细胞测序数据分析
#课程-F:\02.生信课程\01.单细胞2024-可以重点学习
#前部分因为设计多样本导入，采用的是课程
#F:\02.生信课程\01.单细胞2024-可以重点学习\03.单细胞+免疫原性死亡（2024年还需要重点看）\
#04生信培训：同时学会TCGA和单细胞，并整合实战，生信高分必备技能

####清除环境
rm(list=ls())

options(stringsAsFactors = F) #R将字符向量视为字符类型而不是因子类型
#加载R包
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(harmony)
#读取数据
#第一节：批量读取单细胞的数据
#本次采用蜕膜样品进行分析
dir_name=c('GSM6613036',	'GSM6613037',	'GSM6613038',	'GSM6613039',	'GSM6613040',  #正常绒毛
           'GSM6613044',	'GSM6613045',	'GSM6613046' )   #RSA绒毛                   


datalist=list()
for (i in 1:length(dir_name)){
  dir.10x = paste0("GSE214607_RAW/",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  datalist[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i], 
                                   #每个基因至少在3个细胞中表达，每一个细胞至少有250个基因表达
                                   min.cells = 3, min.features = 250)
}
#修改名称
names(datalist)=dir_name


#第二节：细胞质控####
# 批量计算线粒体和rRNA占比
for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 计算线粒体占比
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 计算rRNA占比
  datalist[[i]] <- sce
  rm(sce)
}
#质控前的
violin=list()
for (i in 1:length(datalist)){
  violin[[i]] <- VlnPlot(datalist[[i]],
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                         pt.size = 0.1, 
                         ncol = 4)
}
pearplot_befor <- CombinePlots(plots = violin , nrow=length(datalist), legend="none")
pearplot_befor
#这里图片展示有问题，舍弃
#ggsave(filename = '2.figures/01.QC_before.pdf',plot = pearplot_befor,he=15,wi=15)

#样本合并
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])

#统计每一个样本的个数
raw_count <- table(sce@meta.data$orig.ident)
table(sce@meta.data$orig.ident)

pearplot_befor1<-VlnPlot(sce,
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                         pt.size = 0.1, 
                         ncol = 4)
pearplot_befor1
ggsave(filename = '2.figures/01.QC_before1.pdf',plot = pearplot_befor1,he=10,wi=20)
rm(sce)

#过滤
datalist <- lapply(X = datalist, FUN = function(x) {
  x<-subset(x,subset = nFeature_RNA > 300 & 
              nFeature_RNA < 6000 &  #筛选出每个细胞中检测到的基因数量在500到6000之间的细胞。
              quantile(percent.mt, 0.98) > percent.mt & percent.mt < 25 &
              quantile(percent.Ribo, 0.97) > percent.Ribo & percent.Ribo > quantile(percent.Ribo, 0.01) & 
              nCount_RNA < quantile(nCount_RNA, 0.97) & nCount_RNA > 1000 )
})
#合并数据
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
clean_count <- table(sce@meta.data$orig.ident)
table(sce@meta.data$orig.ident)

#过滤前后样本细胞数据的统计
summary_cells <- as.data.frame(cbind(raw_count,clean_count))
counts <- rbind(as.data.frame(cbind(summary_cells[,1],rep("raw",each = length(summary_cells[,1])))),
                as.data.frame(cbind(summary_cells[,2],rep("clean",each = length(summary_cells[,2])))))
counts$sample <- rep(rownames(summary_cells),times =2)
colnames(counts)<- c("count","Stat","sample")
counts[,1] <- as.numeric(counts[,1])
counts$Stat <- factor(counts$Stat, levels=c("raw", "clean"), ordered=TRUE)
fit_cell_count <- ggplot(data =counts, mapping = aes(x = sample, y=count))+ 
  geom_bar(aes(fill = Stat),stat = 'identity', position = 'dodge') + scale_fill_brewer(palette = "Set4") +
  theme(text=element_text(size=10),legend.title=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black",size = 0.2),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.background = (element_rect(colour= "white",fill = "white")))

fit_cell_count
ggsave(filename = '2.figures/02.fit_cell_count.pdf',plot = fit_cell_count,width = 8,height = 8)

#质控后的小提琴图
violin_after=list()
for (i in 1:length(datalist)){
  violin_after[[i]] <- VlnPlot(datalist[[i]],
                               features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"), 
                               pt.size = 0.1,
                               ncol = 4)
}
pearplot_after <- CombinePlots(plots = violin_after , nrow=length(datalist), legend="none")
pearplot_after
#显示不清，舍弃
#ggsave(filename = '2.figures/04.QC_after.pdf',plot = pearplot_after,he=10,wi=20)

pearplot_after1 <- VlnPlot(sce,
                           features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"), 
                           pt.size = 0.1,
                           ncol = 4)
pearplot_after1
ggsave(filename = '2.figures/03.QC_after1.pdf',plot = pearplot_after1,he=10,wi=20,limitsize = FALSE)
#质控前后图片的合并
pearplot_befor1
pearplot_after1
qc_merge<- CombinePlots(plots = list(pearplot_befor1,pearplot_after1) , 
                        nrow=2, legend='none')
qc_merge
ggsave(filename = '2.figures/04.QC_merge.pdf',plot = qc_merge,he=12,wi=20,limitsize = FALSE)

###从这里开始#课程对接课程：
###   F:\02.生信课程\01.单细胞2024-可以重点学习

#质控是为了去除无效的细胞，质控后再进行合并
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
rm(datalist)


#这里是一步法就做完了，为了导出相关图片，采用下面的分步骤进行
#scRNA_harmony <- NormalizeData(sce) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)

#使用LogNormalize对数据进行归一化
scRNA_harmony <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
rm(sce)
#样本的分组
meta1<-data.frame(matrix(nrow=length(scRNA_harmony@meta.data$orig.ident), ncol=2)) 
#确定行名和顺序
colnames(meta1)=c('Sample','Group1')
meta1$Sample=scRNA_harmony@meta.data$orig.ident
unique(meta1$Sample)
### GSM36-40为 正常绒毛
meta1[grep("GSM6613036",meta1$Sample),]$Group1="normal"
meta1[grep("GSM6613037",meta1$Sample),]$Group1="normal"
meta1[grep("GSM6613038",meta1$Sample),]$Group1="normal"
meta1[grep("GSM6613039",meta1$Sample),]$Group1="normal"
meta1[grep("GSM6613040",meta1$Sample),]$Group1="normal"
### 44-46为 RSA绒毛
meta1[grep("GSM6613044",meta1$Sample),]$Group1="RSA"
meta1[grep("GSM6613045",meta1$Sample),]$Group1="RSA"
meta1[grep("GSM6613046",meta1$Sample),]$Group1="RSA"
#可以导出分组信息
write.table(meta1,file="1.data/1.分组信息.xls",sep="\t",row.names=F,quote=F)

#使用AddMetaData添加分组信息
scRNA_harmony <- AddMetaData(scRNA_harmony, meta1$Sample,col.name = "Sample")
scRNA_harmony <- AddMetaData(scRNA_harmony, meta1$Group1,col.name = "Group1")



##找到高变基因
###官方推荐是2000个高变基因，很多文章也有设置30000的，这个因自己的实验项目决定
scRNA_harmony  <- FindVariableFeatures(scRNA_harmony , selection.method = "vst", nfeatures = 2000) 


#把top20的高变基因挑选出来，目的是为了作图
top20 <- head(VariableFeatures(scRNA_harmony ), 20) 
#画出来不带标签的高变基因图
plot1 <- VariableFeaturePlot(scRNA_harmony ) 
###把top10的基因加到图中
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE, size=2.5,legend="bottom") 
plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom") 
###画图
plot 
ggsave(filename = '2.figures/05.top_20.pdf',plot = plot2,he=8,wi=8)

##如果内存不够，可以只对高变基因进行标准化
#scale.genes <-  VariableFeatures(scRNA)
#scRNA <- ScaleData(scRNA, features = scale.genes)

#对数据进行标准化，占内存
scale.genes <-  rownames(scRNA_harmony)
scRNA_harmony  <- ScaleData(scRNA_harmony , features = scale.genes)


#PCA降维并提取主成分
scRNA_harmony <- RunPCA(scRNA_harmony, features = VariableFeatures(scRNA_harmony),verbose=T) 
plot1 <- DimPlot(scRNA_harmony, reduction = "pca", group.by="orig.ident") 
###画图
plot1 
ggsave(filename = '2.figures/06.sc_pca.pdf',plot = plot1 ,he=8,wi=12)

system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")})

###########这个小步骤就比较关键了，我们需要选择出主成分的数目，用于后续细胞分类。
##值得注意的是，这里定义的“维度”并不代表细胞类型的数目，而是对细胞分类时需要用到的一个参数。
####确定数据的维度 Determine the ‘dimensionality’ of the dataset 
###ElbowPlot() 可以快速的检查降维的效果
plot2 <- ElbowPlot(scRNA_harmony, ndims=50, reduction="pca") 
plot2
#下面的图很浪费时间，所以一般只画上面的图
#scRNA2 <- JackStraw(scRNA_harmony, num.replicate = 100)
#scRNA2 <- ScoreJackStraw(scRNA2, dims = 1:20)
#scRNA2<-JackStrawPlot(scRNA2, dims = 1:20)
#scRNA2
###我们一般选择拐点作为降维的度数。

ggsave("2.figures/07.pca.pdf", plot = plot2, width = 6, height = 6) 
#后续分析要根据右图选择提取的pc轴数量，一般选择斜率平滑的点之前的所有pc轴，此图我的建议是选择前13个pc轴。
##可以看出大概在PC为13的时候，每个轴是有区分意义的。
#自己的数据选择20
pc.num=1:20


#细胞聚类
###一定要指定harmony###这个分辨率是可以自定义的，当我们的样本细胞数较大时候resolution 要高一些，一般情况2万细胞以上都是大于1.0的  
  
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.5)
## 查看每一类有多少个细胞
table(scRNA_harmony@meta.data$seurat_clusters)

length(scRNA_harmony$orig.ident)

##系统发育分析
scRNA_harmony<-BuildClusterTree(scRNA_harmony)
PlotClusterTree(scRNA_harmony)

#pc.num在前面已经指定
#UMAP降维
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = pc.num)
##保存数据这节课的数据
#save(scRNA_harmony,file = '1.data/scRNA_harmony_UMAP降维后.RData')


plot1 =DimPlot(scRNA_harmony, reduction = "umap",label = T) 
plot2 = DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident') 
#combinate
plotc <- plot1+plot2
plotc

ggsave('2.figures/08.sc_umap_按样本和聚类簇展示细胞群.pdf',plotc,he=9,wi=18)

#可视化
sc_umap = DimPlot(scRNA_harmony,
                  #group.by = 'orig.ident',
                  split.by = 'Group1',#决定画图方式
                  reduction="umap",
                  #reduction="tsne",
                  label = "T", 
                  pt.size = 0.2,
                  label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  ) 
sc_umap
ggsave('2.figures/09.sc_umap_按分组展示细胞群.pdf',sc_umap,he=9,wi=18)


##展示自己想要的基因
plota = VlnPlot(scRNA_harmony, 
                features = c("TNRC6B","SRSF3","TBX2"),
                #split.by = 'Group1',
                pt.size = 0,ncol=3)

ggsave('2.figures/10.未注释前展示细胞表达.pdf',plota ,he=9,wi=18)

####单细胞转录组基础分析四：细胞类型鉴定 ####

####细胞类型的注释一般有三种方法
#1、利用marker基因查找网站进行注释  
#2、使用singler进行注释 
#3、根据已有的生物学知识或者文献，按照dotplot来注释。
##现在使用方法一寻找marker基因使用网站注释 找marker基因有以下方法三选一，建议第一种
#默认wilcox方法



#这行代码的意义是将 scRNA_harmony 中的多个数据层（也称为数据层或数据矩阵）合并为一个
#单一的数据层。
# 在单细胞RNA测序数据中，通常会有多个数据层，每个数据层代表不同的生物学特征（例如基因表达）。
# 合并数据层可以将这些信息整合在一起，方便后续的数据分析和处理。
scRNA_harmony <- JoinLayers(scRNA_harmony)



markers <- FindAllMarkers(object = scRNA_harmony, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)   
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top5 = all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5=unique(top5$gene)
sc_marker_dotplot <- DotPlot(object = scRNA_harmony, 
                             features = top5,
                             cols=c("blue", "red"),
                             scale = T)+ 
  RotatedAxis()+ ggtitle("Top 5 Marker Genes")+ 
  theme(plot.title = element_text(hjust = 0.5)) 

sc_marker_dotplot

ggsave(filename = '2.figures/10.sc_marker_dotplot.pdf',
       plot = sc_marker_dotplot,
       height = 9,width = 25)

#热图展示
library(viridisLite)
sc_marker_heatmap<- DoHeatmap(object = scRNA_harmony,
                              features = top5,
                              #group.colors = mycolor,
                              label = F) + 
  ggtitle("Top 5 Marker Genes") + 
  theme(plot.title = element_text(hjust = 0.5)) 
sc_marker_heatmap
ggsave(filename = '2.figures/sc_marker_heatmap.pdf',
       plot = sc_marker_heatmap,
       width = 12,height = 12)



library(SingleR)
load('1.data/hpca.se.Rdata')

ref <- get(load("1.data/BlueprintEncode_bpe.se_human.RData"))
load('1.data/scRNA_harmony.RData')

#获取基因的表达谱的count数据
testdata <- GetAssayData(scRNA_harmony, slot="data")
#获取聚类的亚群
####这里以后可以使用labels = ref$label.main,使用hpca.se就不会出那么多东西


clusters <- scRNA_harmony@meta.data$seurat_clusters
pred.sce <- SingleR(test =  testdata, 
                    ref = ref, 
                    labels = ref$label.fine,
                    method = "clusters",
                    clusters = clusters, 
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")

plotScoreHeatmap(pred.sce)

celltype = data.frame(ClusterID=rownames(pred.sce), celltype=pred.sce$labels
                      , stringsAsFactors = F)
celltype
write.table(celltype,'1.data/celltype.txt',quote = F,sep = '\t',row.names = F)


#修改亚群的名称  根据singleR
scRNA_harmony <- RenameIdents(object = scRNA_harmony, 
                    "0" = celltype[1,2],
                    "1" = celltype[2,2],
                    "2" = celltype[3,2],
                    "3" = celltype[4,2],
                    "4" = celltype[5,2],
                    "5" = celltype[6,2],
                    "6" = celltype[7,2],
                    "7" = celltype[8,2],
                    "8" = celltype[9,2],
                    "9" = celltype[10,2],
                    "10" = celltype[11,2],
                    "11" = celltype[12,2],
                    "12" = celltype[13,2],
                    "13" = celltype[14,2],
                    "14" = celltype[15,2],
                    "15" = celltype[16,2],
                    "16" = celltype[17,2],
                    "17" = celltype[18,2],
                    "18" = celltype[19,2],
                    "19" = celltype[20,2],
                    "20" = celltype[21,2],
                    "21" = celltype[22,2],
                    "22" = celltype[23,2]
                    
)

p222=DimPlot(scRNA_harmony, 
             reduction = "umap",
             pt.size = 0.8,
            label = T,
            split.by = 'Group1',
            label.box = T)
ggsave('2.figures/12.标注细胞群.pdf',p222,he=9,wi=18)


DotPlot(scRNA_harmony,  c("TNRC6B","SRSF3","TBX2"))

#downsample = 100指的是随机抽取100个细胞
DoHeatmap(subset(sce, downsample = 100), features = c("TNRC6B","SRSF3","TBX2"), size = 3)




##展示自己想要的基因
#这里是注释后，因为表达太少几乎看不到，所以采用注释前的图
plota = VlnPlot(scRNA_harmony, 
                features = c("TNRC6B","SRSF3","TBX2"),
                #split.by = 'Group1',
                pt.size = 0,ncol=3)

ggsave('2.figures/998.注释后展示细胞表达.pdf',plota ,he=9,wi=15)



plotb=FeaturePlot(scRNA_harmony,
                  c("TNRC6B","SRSF3","TBX2"),
                  #split.by = 'Group1',
                  pt.size = 0,
                  cols=c("grey",'red'),ncol=3
)

ggsave('2.figures/11.3个基因的FeaturePlot.pdf',plotb ,he=9,wi=27)


plotc=DotPlot(scRNA_harmony,  c("TNRC6B","SRSF3","TBX2"),cols=c("grey",'red'))
ggsave('2.figures/13.3个基因的DotPlot.pdf',plotc ,he=9,wi=9)


#懿芸师姐的分析到此结束


#峰峦图
RidgePlot(scRNA_harmony, 
          features = c("TNRC6B","SRSF3","TBX2"))

#绘制细胞比例柱形图
# 计算scedata数据框中Group1列的频数分布
table(scedata$Group1)

# 计算scedata中每个身份（Idents）的比例
prop.table(table(Idents(scedata)))

# 计算scedata中各组(Group1)不同细胞群(Idents)的细胞数
table(Idents(scedata), scedata$Group1) # 各组不同细胞群的细胞数

# 计算各组样本中不同细胞群的比例，并按列（margin = 2）计算
Cellratio <- prop.table(table(scedata$Group1, Idents(scedata)), margin = 2)
Cellratio
# 将比例表转换为数据框
Cellratio <- as.data.frame(Cellratio)
# 修改Cellratio数据框的第一和第二列列名分别为Group和Celltype
colnames(Cellratio)[1:2] <- c("Group", "Celltype")
Cellratio
# 计算不同组(Group)的数量，用于设置颜色数量
colourCount = length(unique(Cellratio$Group))
# 加载ggplot2包
library(ggplot2)
# 创建一个ggplot对象，以Cellratio数据框为数据源
plotd = ggplot(Cellratio) + 
  # 使用geom_bar()函数绘制条形图
  # aes()函数设置美学映射，x轴为Celltype，y轴为Freq，填充颜色根据Group分组
  geom_bar(aes(x = Celltype, y = Freq, fill = Group), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') + 
  # 使用经典主题
  theme_classic() +
  # 设置x轴和y轴标签
  labs(x = 'Celltype', y = 'Ratio') +
  # 翻转坐标轴
  coord_flip() +
  # 设置面板边框的样式
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        # 修改Y轴标签的字体大小
        axis.text.y = element_text(size = 14),
        # 修改X轴坐标值的字体大小
        axis.text.x = element_text(size = 14),
        # 修改X轴和Y轴标题的字体大小
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        # 修改图例的位置到上方
        legend.position = "top",
        # 修改图例字体大小
        legend.text = element_text(size = 14)) +
  # 手动设置填充颜色，从#FF9999到#9999FF
  scale_fill_manual(values = c("#FF9999", "#9999FF"))

# 打印图形
print(plotd)


ggsave('2.figures/14.细胞比例柱状图.pdf',plotd ,he=9,wi=9)

save(scRNA_harmony,file = '1.data/scRNA_harmony.RData')


#一些常用的标记物
genes_to_check = c('EPCAM','KRT19','CLDN4',  #上皮
                   'PECAM1' , 'CLO1A2', 'VWF',  #基质
                   'CD3D', 'CD3E', 'CD8A', 'CD4','CD2', #T
                   'CDH5', 'PECAM1', 'VWF',  #内皮
                   'LUM' , 'FGF7', 'MME',  #成纤维
                   'AIF1', 'C1QC','C1QB','LYZ',  #巨噬
                   'MKI67', 'STMN1', 'PCNA',  #增殖
                   'CPA3' ,'CST3', 'KIT', 'TPSAB1','TPSB2',#肥大
                   'GOS2', 'S100A9','S100A8','CXCL8', #中性粒细胞
                   'KLRD1', 'GNLY', 'KLRF1','AREG', 'XCL2','HSPA6', #NK
                   'MS4A1','CD19', 'CD79A','IGHG1','MZB1', 'SDC1',  #B
                   'CSF1R', 'CSF3R', 'CD68') #髓系




p333 = DotPlot(scRNA_harmony, features = unique(genes_to_check),
               assay='RNA'  )  + coord_flip()

p333


VlnPlot(scRNA_harmony, 
        features = c('LUM' , 'FGF7', 'MME'),
        #split.by = 'Group1',
        pt.size = 0,ncol=3)
