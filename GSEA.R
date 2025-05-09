#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

allFile="all.txt"      #所有基因的差异结果文件
gmtFile="c2.cp.kegg.symbols.gmt"      #基因集文件
#setwd("C:\\Users\\lexb\\Desktop\\ICD\\20.GSEA")     #设置工作目录

#读取输入文件,并对输入文件进行整理
rt=read.table(allFile, header=T, sep="\t", check.names=F)
rt=rt[order(rt[,"logFC"],decreasing=T),]
logFC=as.vector(rt[,"logFC"])
names(logFC)=as.vector(rt[,1])

#读入基因集文件
gmt=read.gmt(gmtFile)

#对排序好的基因进行GSEA富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)

#输出ICD高表达组富集的图形
kkUp=kkTab[kkTab$NES>0,]
termNum=5     #设置展示通路的数目，展示前5个富集最显著的通路
if(nrow(kkUp)>=termNum){
  showTerm=row.names(kkUp)[1:termNum]      #获取展示通路的名称
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in RSA")
  pdf(file="GSEA.RSA.pdf", width=7, height=5.5)
  print(gseaplot)
  dev.off()
}

#输出ICD低表达组富集的图形
kkDown=kkTab[kkTab$NES<0,]
termNum=5     #设置展示通路的数目，展示前5个富集最显著的通路
if(nrow(kkDown)>=termNum){
  showTerm=row.names(kkDown)[1:termNum]      #获取展示通路的名称
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in Control")
  pdf(file="GSEA.Control.pdf", width=7, height=5.5)
  print(gseaplot)
  dev.off()
}


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: eduBio


