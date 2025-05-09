
#加载 limma 和 sva 文库，其中包含分析和校正基因表达数据的功能。
library(limma)
library(sva)

#这一行定义了一个变量 outFile，并将其设置为“ merge.txt”，
#这是要写入合并数据的文件的名称
outFile="merge.txt"       


#这两行使用 dir ()函数列出工作目录中的所有文件，然后过滤列表，
#只包含带有“。“ txt”扩展名。生成的文件名列表存储在 files 变量中。
files=dir()
files=grep("txt$", files, value=T)

#这些行初始化一个空列表 genList，然后循环遍历文件中的每个文件。
#对于每个文件，脚本使用 read.table ()从文件中读取数据，从数据的第一列中提取基因名称
#使用 only ()删除重复，然后使用文件名前缀(通过将文件名分割到“”上获得)
#将得到的唯一基因名称列表存储在 genList 列表中和“-”)作为列表中的键。
geneList=list()
for(file in files){
	if(file==outFile){next}
    rt=read.table(file, header=T, sep="\t", check.names=F)     
    geneNames=as.vector(rt[,1])      
    uniqGene=unique(geneNames)       
    header=unlist(strsplit(file, "\\.|\\-"))
    geneList[[header[1]]]=uniqGene
}

#这一行使用 Reduce ()函数将 intersect ()函数应用于 genList 列表中的所有元素。
#这导致出现在所有输入文件中的基因名称向量
interGenes=Reduce(intersect, geneList)

#这些行初始化一个空数据帧 allTab 和一个空向量 batchType。
#这些变量将分别用于存储组合的基因表达数据和批处理信息
allTab=data.frame()
batchType=c()

#以表的形式读入文件，指定第一行为标题，分隔符为制表符。
#提取基因名称作为载体。
#查找唯一的基因名称并将它们保存为名为 uniqGene 的载体。
#提取文件名，并使用它为 genList 列表对象创建一个头
#将 genList 对象中的 uniqGene 向量保存在适当的头下。
for(i in 1:length(files)){
    inputFile=files[i]
    header=unlist(strsplit(inputFile, "\\.|\\-"))
    rt=read.table(inputFile, header=T, sep="\t", check.names=F)
    rt=as.matrix(rt)
    rownames(rt)=rt[,1]
    exp=rt[,2:ncol(rt)]
    dimnames=list(rownames(exp),colnames(exp))
    data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
    rt=avereps(data)

    #判断是否需要log2转换
    qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
    if(LogC){
    	rt[rt<0]=0
        rt=log2(rt+1)}
    rt=normalizeBetweenArrays(rt)
    
   
    if(i==1){
    	allTab=rt[interGenes,]
    }else{
    	allTab=cbind(allTab, rt[interGenes,])
    }
    batchType=c(batchType, rep(i,ncol(rt)))
}

#最后一步使用 sva 包中的 ComBat 函数对 allTab 数据执行批量校正。
#BatchType 变量用于标识每个示例的批号。
#然后，将生成的数据帧写入名为 merge.txt 的文件。
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge.txt", sep="\t", quote=F, col.names=F)


