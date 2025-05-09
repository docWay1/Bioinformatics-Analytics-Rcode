# 导入所需的包
library(xgboost)
library(data.table)
library(ggplot2)

set.seed(123)
#读取输入文件
data=read.table("vennGeneExp.txt" , header=T, sep="\t", check.names=F, row.names=1)
# 创建一个数据框，用于构建xgboost模型
group = as.factor(data$group)
# 提取特征数据（除去 'group' 列）
data = data[, -which(names(data) %in% c("group"))]


df <- data.frame(data)
# 提取分组变量并转化为因子






# 准备训练数据
X <- df
y <- group

# 将数据转换为DMatrix格式，以便用于xgboost模型
dtrain <- xgb.DMatrix(data = as.matrix(X), label = as.integer(y) - 1)  # xgboost使用类别从0开始编号

# 定义xgboost参数
params <- list(
  objective = "multi:softprob",  # 多分类问题的目标函数
  num_class = length(levels(y)),  # 类别数量
  eval_metric = "mlogloss"  # 使用多分类对数损失进行模型评估
)

# 训练xgboost模型
model <- xgb.train(params = params, data = dtrain, nrounds = 1000)

# 获取特征重要性
feature_importance <- xgb.importance(model = model)
print(feature_importance)

# 选择重要性较高的特征基因
top_features <- feature_importance[order(-feature_importance$Cover), ]

# 选择重要性较高的特征基因（根据Gain>0.01）
top_features <- feature_importance[feature_importance$Gain > 0.00, ]

# 输出结果
print(top_features)


write.table(top_features$Feature, file="xgboost.gene_selected.txt", sep="\t", quote=F, row.names=F, col.names=F)

write.table(feature_importance, file="xgboost.result.txt", sep="\t", quote=F, row.names=F, col.names=T)




# 特征重要性图（在X=0.01处添加一条黑色虚线和文本标签）
p1 <- ggplot(data = feature_importance, aes(y = reorder(Feature, Gain), x = Gain, fill = Gain)) +
  geom_bar(stat = "identity") +
  #geom_vline(xintercept = 0.01, linetype = "dashed", color = "black") +  # 添加虚线
  geom_text(aes(label = "0.3"), x = 0.3, y = Inf, vjust = -0.5, hjust = 0, color = "black") +  # 添加文本标签
  labs(x = "Importance", y = "Feature",title = "Feature gene importance (XGBoost)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank()) +
  scale_fill_gradient(low = "pink", high = "red")




pdf(file="xgboost_Gain.pdf", width=6, height=6)
p1
dev.off()
p1
