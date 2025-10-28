# AutoExon

### 更新日志

#### v2.0.0

- date: 2025/10/28
- 新结构

#### v1.2.2

- date: 2025/10/28
- 新的psm计算逻辑，以pepmapviz坐标为标准

#### v1.2.2

- date: 2025/10/28
- 新的psm计算逻辑，以uniport坐标为标准

#### v1.2.1

- date: 2025/10/21
- 修改保存文件逻辑

#### v1.2.0

- date: 2025/10/01
- 修复数据映射BUG
- 修复坐标overlap检测BUG
- 修复排序BUG

#### v1.1.0

- date: 2025/10/01
- 加入多进程
- TTN处理有问题

#### v1.0.0

- date: 2025/10/01
- 流程正常
- TTN处理有问题

### R环境安装

```
apt-get update && apt-get install -y r-base r-base-dev
Rscript -e 'install.packages(c("ggplot2"), repos="https://cloud.r-project.org")'
Rscript -e 'install.packages(c("data.table","dplyr","jsonlite","magrittr","remotes"), repos="https://cloud.r-project.org")'
Rscript -e 'remotes::install_github("Genentech/PepMapViz")'
```

