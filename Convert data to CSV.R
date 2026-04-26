# 如果没装过，先装（只需一次）
# install.packages(c("jsonlite", "data.table"))

library(jsonlite)
library(data.table)

# 1. 设置工作目录（注意：R 里用 / 或 \\）
setwd("D:/桌面临时文件/已经发表论文/慧玲论文/FDA微生物/分析数据")

# 2. 读取 JSON
json_data <- fromJSON(
  "food-enforcement-0001-of-0001.json",
  flatten = TRUE
)

# 3. 提取真正的数据部分（results）
df <- as.data.table(json_data$results)

# 4. 导出为 CSV（Excel / Stata 友好）
fwrite(
  df,
  file = "food_enforcement.csv",
  bom = TRUE   # 防止 Excel 中文乱码
)

# 5. 简单检查
dim(df)
names(df)
