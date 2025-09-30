# AutoExon

### 更新日志

#### v1.0

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

