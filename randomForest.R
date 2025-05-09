

# 引用包
library(randomForest)
set.seed(127)  
inputFile="vennGeneExp.txt"       # 输入文件

# 读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

# 从数据中提取分组变量并作为响应变量
group = as.factor(data$group)

# 提取除了 'id' 和 'group' 的所有列作为特征变量
data = data[, -which(names(data) %in% c("group"))]

# 随机森林模型
rf=randomForest(group~., data=data, ntree=2500)
pdf(file="randomforest.pdf", width=6, height=6)
plot(rf, main="Random Forest Error Rate", lwd=2)
dev.off()

# 找出误差最小的点
optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(group~., data=data, ntree=optionTrees)

# 查看基因的重要性
importance=importance(x=rf2)

# 绘制基因的重要性图
pdf(file="geneImportance.pdf", width=6, height=6)
varImpPlot(rf2, main="Gene Importance")
dev.off()

# 挑选疾病特征基因
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes

# 挑选重要性评分大于0的基因
rfGenes=names(rfGenes[rfGenes > 0.5])     
write.table(rfGenes, file="rfGenes_selected.txt", sep="\t", quote=F, col.names=F, row.names=F)

