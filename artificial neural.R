# 引入必要的库
library(limma)               # 用于微阵列数据分析
library(neuralnet)           # 用于训练神经网络
library(NeuralNetTools)      # 提供绘制神经网络的工具
library(pROC)                # 用于生成ROC曲线

# 设置种子以确保结果的可重复性
set.seed(123)

# 文件路径设置
expFile="3GeneExp.txt"     
diffFile="diff.txt"        
inputFile="geneScore.txt"

# 读取和处理表达量数据
rt = read.table(expFile, header=TRUE, sep="\t", check.names=FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[,1]
exp = rt[,2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data = avereps(data)  # 平均重复实验的数据

# 读取差异表达数据
diffRT = read.table(diffFile, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
diffRT = diffRT[row.names(data),]

# 分类上调和下调的基因数据
dataUp = data[diffRT["logFC"] > 0,]
dataUp2 = t(apply(dataUp, 1, function(x) ifelse(x > median(x), 1, 0)))

# 保存评分数据
write.table(dataUp2, file=inputFile, sep="\t", quote=FALSE, col.names=NA)

# 读取评分文件并准备数据
data=read.table(inputFile, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
data=as.data.frame(t(data))

# 获取样本分组信息
group=gsub("(.*)_(.*)", "\\2", row.names(data))
data$con=ifelse(group=="con", 1, 0)
data$treat=ifelse(group=="treat", 1, 0)

# 训练神经网络模型
fit=neuralnet(con+treat~., data, hidden=5)

# 绘制神经网络结构并保存为PDF
pdf(file="neuralnet2.pdf", width=9, height=9)
plotnet(fit)
dev.off()

# 利用模型进行预测并生成预测准确率表
net.predict=compute(fit, data)$net.result
net.prediction=c("con", "treat")[apply(net.predict, 1, which.max)]
predict.table=table(group, net.prediction)

# 计算预测准确性
conAccuracy=predict.table[1,1]/(predict.table[1,1]+predict.table[1,2])
treatAccuracy=predict.table[2,2]/(predict.table[2,1]+predict.table[2,2])

# 输出预测结果和准确性
write.table(net.predict, file="neural.predict.txt", sep="\t", quote=FALSE, col.names=FALSE)

# 读取预测结果文件，准备ROC分析数据
rt=read.table("neural.predict.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)
y=ifelse(gsub("(.*)_(.*)", "\\2", row.names(rt))=="con", 0, 1)

# 生成ROC曲线并保存为PDF
roc1=roc(y, as.numeric(rt[,2]))
ci1=ci.auc(roc1, method="bootstrap")
ciVec=as.numeric(ci1)


pdf(file="ROC.pdf", width=6, height=6)
#plot(roc1, print.auc=TRUE, col="red", legacy.axes=TRUE, main="")
#text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f",ciVec[1]), "-", sprintf("%.03f",ciVec[3])), col="red")
# 绘制ROC曲线，颜色设为红色，不自动打印AUC值，使用传统坐标轴样式，主标题为基因名
plot(roc1, col="skyblue", legacy.axes=T, main="ANN model")
# 手动添加AUC值的文本，颜色设为黑色
text(0.23, 0.06, paste0("AUC: ", sprintf("%.03f", roc1$auc)), col="red")
# 添加95% CI的文本，颜色设为黑色
text(0.15, 0.01, paste0("95% CI: ", sprintf("%.03f", ciVec[1]), "-", sprintf("%.03f", ciVec[3])), col="red")
dev.off()

# 打印准确率
print(paste0("Con accuracy: ", sprintf("%.3f", conAccuracy)))
print(paste0("Treat accuracy: ", sprintf("%.3f", treatAccuracy)))
