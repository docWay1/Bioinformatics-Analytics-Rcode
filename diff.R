
library(limma)
library(pheatmap)

inputFile="merge.txt"       
logFCfilter=0.3        
adj.P.Val.Filter=0.05       



rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#这部分代码从工作目录中分别以“s1.txt”和“s2.txt”结尾的两个文件中读取示例的名称。
#然后，它将每个文件中的唯一基因名称组合到两个单独的载体中，“sampleName1”和“sampleName2”。
sampleName1=c()
files=dir()
files=grep("s1.txt$", files, value=T)
for(file in files){
    rt=read.table(file, header=F, sep="\t", check.names=F)      
    geneNames=as.vector(rt[,1])      
    uniqGene=unique(geneNames)       
    sampleName1=c(sampleName1, uniqGene)
}


sampleName2=c()
files=dir()
files=grep("s2.txt$", files, value=T)
for(file in files){
    rt=read.table(file, header=F, sep="\t", check.names=F)      
    geneNames=as.vector(rt[,1])     
    uniqGene=unique(geneNames)       
    sampleName2=c(sampleName2, uniqGene)
}

#这些代码行根据分别以“s1.txt”和“s2.txt”结尾的文件中提供的样本名称，
#将输入数据矩阵数据拆分为两个单独的矩阵 conData 和 treatData。

#然后使用 cbind 函数将两个矩阵 conData 和 treatData 逐列连接，
#以形成具有来自对照组和处理组的样本的新矩阵数据。conData 和 treatData 的列数
#也分别存储为变量 conNum 和 treatNum。
conData=data[,sampleName1]
treatData=data[,sampleName2]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

#这行代码是在创建一个类型向量(Type vector)，
#其中向量中的元素数量等于对照组(con)样本的数量(conNum)加上治疗组(treat)样本的数量(treatNum)。
#其中，向量中前conNum个元素为"con"，后treatNum个元素为"treat"。
#这个向量将在后面的代码中用于定义线性模型的设计矩阵。
Type=c(rep("con",conNum),rep("treat",treatNum))
#首先，代码将数据集中两个样本的数据提取出来，并用 cbind 函数将它们合并成一个数据集。然后根据样本的类别（对照组或实验组）生成一个设计矩阵 design。
#接着，使用线性模型 lmFit 对数据进行拟合，并使用 makeContrasts 函数生成对照矩阵 cont.matrix，用于比较实验组和对照组之间的差异。
#最后，使用 eBayes 函数对差异分析进行贝叶斯校正，得到基因的表达差异结果。
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#adjust='fdr'表示使用FDR校正方法，number=200000表示返回前200000个基因。
#topTable函数将根据提供的参数对差异进行排序，并将结果存储在allDiff中。
#最后，将基因名称添加到allDiff的顶部，以便在将结果写入文件时可以包括基因名称。
allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)

#这段代码的作用是将数据和样本类型合并
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)

#选择使用P.Value来进行分析
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & P.Value < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)




diffGeneExp=data[row.names(diffSig),]
diffGeneExpOut=rbind(id=paste0(colnames(diffGeneExp),"_",Type), diffGeneExp)
write.table(diffGeneExpOut, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)




geneNum=50
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
hmExp=data[hmGene,]
Type=c(rep("Con",conNum),rep("Treat",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=10, height=8)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("#9999FF", "white", "#FF9999"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=7,
         fontsize_col=8)
dev.off()


####火山图

library(dplyr)
library(ggplot2)
library(ggrepel)

rtVol = read.table("all.txt", header=T, sep="\t", check.names=F)

Sig=ifelse((rtVol$P.Value<adj.P.Val.Filter) & (abs(rtVol$logFC)>logFCfilter), ifelse(rtVol$logFC>logFCfilter,"Up","Down"), "Not")


rtVol = mutate(rtVol, Sig=Sig)
p = ggplot(rtVol, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(col=Sig))+
  scale_color_manual(values=c("green", "black","red"))+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))

p1=p+geom_label_repel(data=filter(rtVol, ((rtVol$P.Value<adj.P.Val.Filter) & (abs(rtVol$logFC)>logFCfilter))),
                      box.padding=0.1, point.padding=0.1, min.segment.length=0.05,
                      size=1.8, aes(label=id)) + theme_bw()

pdf(file="vol.pdf", width=7, height=6.1)
print(p1)
dev.off()

