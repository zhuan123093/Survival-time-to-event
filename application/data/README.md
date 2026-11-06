# ACTG175 数据来源

`ACTG175.txt` 是 R 包 `speff2trial` 提供的2139人、27变量数据，原文件直接保留。`ACTG175.Rd` 为同一版本的变量和终点说明。来源固定在 CRAN 镜像提交 `9fa954e262daa0dd8d078d8a843537467f7603df`：

- [数据文件](https://github.com/cran/speff2trial/blob/9fa954e262daa0dd8d078d8a843537467f7603df/data/ACTG175.txt)
- [变量说明](https://github.com/cran/speff2trial/blob/9fa954e262daa0dd8d078d8a843537467f7603df/man/ACTG175.Rd)

`sources.json` 记录各原文件的下载地址、大小和SHA256。分析脚本从完整数据选取 Z3 中的两个单药治疗臂；入选行及受试者编号写入 `application/results/sample_index.csv`。

临床试验出处：Hammer et al. (1996), *New England Journal of Medicine*, 335(15):1081–1090, [doi:10.1056/NEJM199610103351501](https://doi.org/10.1056/NEJM199610103351501)。
