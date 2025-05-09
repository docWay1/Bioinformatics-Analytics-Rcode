# 安装Boruta包，如果未安装可以取消注释
#install.packages("Boruta")
library(Boruta) # 加载Boruta包
set.seed(123) # 设置随机数种子，确保结果的可重复性

# 读取数据
BreastCancer=read.table("vennGeneExp.txt", header=T, sep="\t", check.names=F)
str(BreastCancer) # 显示数据结构

# 数据预处理
data = na.omit(BreastCancer) # 删除含有缺失值的行
data = data[, -1] # 删除第一列（如果第一列是无关的标识列）
sum(is.na(data))  # 检测数据是否还有缺失值
data$group = factor(data$group) # 将诊断结果转化为因子类型，适用于分类

# 使用Boruta算法进行特征选择
Boruta.srx <- Boruta(group ~ ., data = data, maxRuns=1000)
print(Boruta.srx) # 打印Boruta算法的结果

# 获取选择的属性
getSelectedAttributes(Boruta.srx, withTentative = F) # 获取最终确认的特征
# 获取包含已确认特征的模型公式
getConfirmedFormula(Boruta.srx) # 获取基于确认特征的模型公式
#导出筛选的基因
write.table(file="Boruta.gene_selected.txt", getSelectedAttributes(Boruta.srx, withTentative = F) , sep="\t", quote=F, row.names=F, col.names=F)

# 获取特征的统计数据
attStats(Boruta.srx) # 获取特征的统计分析结果

# 生成PDF文件，记录特征的重要性得分
pdf("01.各个变量的重要性得分.pdf", width = 6, height = 6)
plot(Boruta.srx, las = 2,xlab = "",cex.axis = 0.7,cex= 0.3) # 绘制特征重要性的图表
dev.off() # 关闭PDF设备

# 分析是否有必要增加迭代次数
# 生成PDF文件，查看Boruta算法运行期间属性Z-Score的变化
pdf("02.Boruta运行期间属性Z-Score的演变.pdf", width = 6, height = 6)
plotImpHistory(Boruta.srx) # 绘制属性Z-Score的历史变化图
dev.off() # 关闭PDF设备

