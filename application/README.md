# ACTG175 Z3 实际应用

分析既往抗病毒治疗超过52周的受试者，比较 didanosine（A=1）与 zidovudine（A=0）。筛选规则为 `strat == 3 & arms %in% c(0, 3)`，共434人、148事件；两组人数分别为221和213，事件分别为56和92。

终点为首次出现CD4计数下降至少50%、进展为AIDS的事件或死亡，观察时间和事件指示分别为 `days` 和 `cens`。数据文件和27个变量的原始说明见 [data/](data/README.md)。

## 调整变量

p不含截距，三组基线变量依次嵌套：

|p|变量|
|---|---|
|5|`cd40`, `preanti`, `cd80`, `age`, `wtkg`|
|8|再加入 `karnof`, `symptom`, `gender`|
|14|再加入 `hemo`, `homo`, `drugs`, `race`, `oprior`, `z30`|

所有设置使用同一批434人。连续变量以线性项进入模型，二元变量使用原始0/1编码；协变量整体中心化并标准化。

## 复算

在仓库根目录运行，R依赖 `survival` 和 `MASS`：

```powershell
Rscript --vanilla application/code/run.R
```

脚本调用根目录 `code/data.R`、`code/methods.R` 和 `code/chiang_loo.R`，计算 lr、ols、dof、loo、cloo、ploo 六种方法。`ols` 对应YE调整检验；`loo` 使用组内含截距CF分子和Chiang mHC3。并列事件在同一事件时点合并，lr 使用 `survival::survdiff` 的并列事件方差修正，调整项使用实际治疗比例221/434。

结果写入 `application/results/`，也可通过 `--output=目录` 指定输出位置。

|文件|内容|
|---|---|
|`test_statistics.csv`|18个方法设置的精确分子、方差、标准差、统计量，以及PLOO/CLOO的I、II|
|`table.csv`, `table.tex`|论文应用表的数值及LaTeX表格；TeX需加载booktabs|
|`sample_index.csv`, `sample_summary.csv`|入选数据行、受试者编号及治疗组人数/事件|
|`design_diagnostics.csv`|变量列表、组内设计秩及最大杠杆值|
|`formula_checks.csv`|派生响应得分恒等式及PLOO/CLOO密集公式核对|

表中分子为 sqrt(n)U，标准差为sigma，统计量为sqrt(n)U/sigma。

原始试验报告为 [Hammer et al. (1996)](https://doi.org/10.1056/NEJM199610103351501)，BibTeX条目见 [reference.bib](reference.bib)。
