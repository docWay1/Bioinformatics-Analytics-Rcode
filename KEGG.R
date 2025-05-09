######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")
#install.packages("digest")
#install.packages("GOplot")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GOplot)

pvalueFilter=0.05       #p值过滤条件
qvalueFilter=0.05       #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

#setwd("C:\\biowolf\\neuralDiagnostic\\10.KEGG")      #设置工作目录
rt=read.table("差异表达CSG.txt", header=T, sep="\t", check.names=F)      #读取输入文件

#基因名字转换为基因id
colnames(rt)[1]="Gene"
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$Gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存显著富集的结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#定义显示通路的数目
showNum=20
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

#柱状图
pdf(file="barplot.pdf", width=8, height=7)
#barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
dev.off()

#气泡图
pdf(file="bubble.pdf", width=8, height=7)
#dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
dev.off()

#获取通路信息
kegg=data.frame(Category="ALL", ID = KEGG$ID, Term=KEGG$Description, Genes = gsub("/", ", ", KEGG$geneID), adj_pval = KEGG$p.adjust)
#读取基因的差异信息
genelist <- data.frame(ID=rt$Gene, logFC=rt$logFC)
row.names(genelist)=genelist[,1]
#设置圈图参数
circ <- circle_dat(kegg, genelist)
termNum =8       #显示通路的数目
termNum=ifelse(nrow(kegg)<termNum,nrow(kegg),termNum)
geneNum=200      #显示基因的数目
geneNum=ifelse(nrow(genelist)<geneNum, nrow(genelist), geneNum)
#绘制通路的圈图
chord <- chord_dat(circ, genelist[1:geneNum,], kegg$Term[1:termNum])
pdf(file="KEGGcircos.pdf", width=10, height=10)
GOChord(chord, 
        space = 0.001,           #基因之间的间距
        gene.order = 'logFC',    #按照logFC值对基因排序
        gene.space = 0.25,       #基因名跟圆圈的相对距离
        gene.size = 5,           #基因名字体大小 
        border.size = 0.1,       #线条粗细
        process.label = 6)       #通路字体大小
dev.off()

#通路的聚类图
pdf(file="KEGGcluster.pdf",width=12, height=10)
GOCluster(circ, 
          kegg$Term[1:termNum], 
          lfc.space = 0.2,        #logFC与树之间的空隙大小
          lfc.width = 1,          #logFC的圆圈宽度
          term.space = 0.2,       #logFC与通路间的空隙大小
          term.width = 1)         #通路圆圈的宽度
dev.off()          




