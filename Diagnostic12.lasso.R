# 设置随机种子
set.seed(123)
library(glmnet)                    

# 读取数据文件
inputFile="vennGeneExp.txt"       

# 读取表格数据，指定文件、表头、分隔符、行名
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

# 提取分组变量作为响应变量
y = rt$group

# 提取表达数据（除去id和group列）
x = rt[, -which(names(rt) %in% c("group"))]

# 将表达数据转换为矩阵
x=as.matrix(x)

# 进行LASSO回归分析
fit=glmnet(x, y, family = "binomial", alpha=1)
cvfit=cv.glmnet(x, y, family="binomial", alpha=1, type.measure='deviance', nfolds = 10)

# 画出交叉验证图
pdf(file="cvfit交叉验证图.pdf",width=6,height=6)
plot(cvfit)
dev.off()

# 画出LASSO系数分布图
pdf(file="lasso系数分布图.pdf",width=6,height=6)
plot(fit, label=TRUE, xlab="L1 Norm", ylab="Coefficients", main="LASSO Coefficient Distribution")
dev.off()

# 筛选并输出非零系数的基因
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]  # 移除第一个元素，通常是截距
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)
lassoGene
