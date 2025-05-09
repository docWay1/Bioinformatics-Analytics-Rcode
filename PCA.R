######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")


#???ð?
library(limma)
library(ggplot2)

clusterFile="PCAdata.txt"     #???͵Ľ????ļ?
#setwd("C:\\biowolf\\m6A\\17.PCA")      #???ù???Ŀ¼

#??ȡ?????ļ?,?????????ļ?????????
rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
data=rt[,1:(ncol(rt)-1),drop=F]
Sample=as.vector(rt[,ncol(rt)])

#PCA????
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Sample=Sample)
PCA.mean=aggregate(PCA[,1:2], list(Sample=PCA$Sample), mean)

#??????ɫ
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
m6aCluCol=bioCol[1:length(levels(factor(Sample)))]


#??????Բ????
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(factor(PCA$Sample))){
df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$Sample==g,],
                  veganCovEllipse(cov.wt(cbind(PC1,PC2),
                  wt=rep(1/length(PC1),length(PC1)))$cov,
                  center=c(mean(PC1),mean(PC2))))), Sample=g))
}

#????PCAͼ??
pdf(file="PCA全部基因.pdf", height=5, width=6.5)
ggplot(data = PCA, aes(PC1, PC2)) + 
  geom_point(aes(color = Sample), size = 4) +  # 设置点的大小为 4（可以根据需要调整）
  scale_colour_manual(name = "Sample", values = m6aCluCol) +
  theme_bw() +
  theme(plot.margin = unit(rep(1.5, 4), 'lines')) +
  geom_path(data = df_ell, aes(x = PC1, y = PC2, colour = Sample), size = 1, linetype = 2) +
  annotate("text", x = PCA.mean$PC1, y = PCA.mean$PC2, label = PCA.mean$Sample, cex = 7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

