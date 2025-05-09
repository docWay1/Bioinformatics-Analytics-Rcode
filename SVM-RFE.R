library(e1071)
library(kernlab)
library(caret)

set.seed(2026)  # 设定随机种子以保证结果的可重复性

inputFile="vennGeneExp.txt"  # 输入文件

# 读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

# 提取分组变量并转化为因子
group = as.factor(data$group)

# 提取特征数据（除去 'group' 列）
data = data[, -which(names(data) %in% c("group"))]

# SVM-RFE 特征选择
Profile = rfe(x=data,
              y=as.numeric(group),
              sizes = c(2,4,6,8, seq(10,40,by=3)),
              rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
              methods="svmRadial")



# 绘制SVM-RFE误差图
pdf(file="SVM-RFE_误差图.pdf", width=6, height=6)
par(las=1)
x = Profile$results$Variables
y = Profile$results$RMSE  # 使用您指定的RMSE作为准确性的反映
plot(x, y, xlab="Number of Features", ylab="RMSE (Cross-Validation)", col="red")
lines(x, y, col="red")
dev.off()

# 绘制SVM-RFE准确性图
pdf(file="SVM-RFE_准确性.pdf", width=6, height=6)
par(las=1)
x = Profile$results$Variables
y = 1 - Profile$results$RMSE  # 误差计算为1-RMSE
plot(x, y, xlab="Number of Features", ylab="1-RMSE (Cross-Validation)", col="red")
lines(x, y, col="red")
dev.off()

# 导出选定的特征基因
featureGenes = Profile$optVariables
write.table(featureGenes, file="SVM-RFE_gene_selected.txt", sep="\t", quote=F, row.names=F, col.names=F)

featureGenes




