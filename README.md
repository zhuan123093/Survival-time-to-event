# 协变量调整 log-rank 检验：模拟与实际应用

本项目比较协变量调整 log-rank 检验在完全随机分配（CRE）、Bernoulli分配（BRE）和分层随机分配（CAR）下的表现，提供七种方法的模拟代码、结果及九张实验图，并提供 ACTG175 Z3 实际应用的数据、六方法分析代码和结果表。

```text
Survival-logrank-reproducibility/
├── README.md                 # 本文件：代码用途、数据说明、运行命令
├── INDEX.csv                 # 九图逐一对应图片、绘图CSV、原始结果
├── MANIFEST.csv              # 文件大小和SHA256
├── code/
│   ├── settings.csv          # 24组设置：n、CRE/BRE/CAR、H0/H1、比例、p、试验偏移
│   ├── data.R                # 70%事件删失、风险集派生响应、标准log-rank
│   ├── methods.R             # OLS/dof、CLOO/PLOO/oracle及共用矩阵运算
│   ├── chiang_loo.R          # 含截距CF分子、Chiang mHC3、饱和组及大样本计算
│   ├── experiment.R          # 资料生成、随机种子、七方法单次计算
│   ├── run.R                 # 按设置运行，或重放已保存记录核对代码
│   ├── summarize.R           # 从24份raw收集汇总，或核对保存的summary
│   └── plot.py               # 绘制九张实验图，导出各图CSV
├── results/
│   ├── raw/                  # 24份RDS，每组10000次模拟的逐次统计量
│   ├── summary.csv           # 700个方法设置的汇总
│   ├── variance_components.csv # LOO方差分项
│   ├── pooled_variance_components.csv # CLOO/PLOO/oracle的I、II和层间项
│   ├── invalid_trials.csv    # 七种方法的无效统计量记录
│   ├── figures/              # 九张图的PDF和PNG
│   ├── plot_data/            # 每图一份CSV，共1134个绘图点
│   ├── code_check.csv        # 重放记录与保存结果的比较
│   ├── summary_check.csv     # raw重计与summary逐列核对
│   └── plot_check.csv        # 绘图数据和图片检查
└── application/
    ├── README.md             # 亚组、终点、调整变量及运行方法
    ├── data/                 # ACTG175完整数据、原始变量说明、来源清单
    ├── code/run.R            # Z3、p=5/8/14、六种检验及表格导出
    ├── results/              # 精确统计量、样本清单、CSV与LaTeX表
    └── reference.bib         # Hammer et al. (1996)原试验词条
```

## 代码对应关系

`run.R` → `experiment.R`：按 `settings.csv` 生成资料，然后调用 `data.R` 计算风险集响应，调用 `methods.R` 和 `chiang_loo.R` 计算七种方法。全部依赖在本目录内。

`summarize.R` → `results/raw/` → `summary.csv`；`plot.py` → `summary.csv` → 九张图片及各图CSV。

|图片文件名（PDF/PNG同名）|内容|
|---|---|
|`introduction_CRE`|引言：CRE七方法比较|
|`CRE_main`|正文：CRE|
|`CAR_main`|正文：CAR|
|`BRE_appendix`|附录：BRE|
|`CRE_with_ols`|附录：CRE含OLS|
|`BRE_with_ols`|附录：BRE含OLS|
|`CAR_with_ols`|附录：CAR含OLS|
|`sample_size_appendix`|附录：固定p=20、增加n|
|`gamma_CRE_appendix`|附录：Gamma协变量CRE|

## 模拟设置和数据

七种方法为lr、ols、dof、loo、cloo、ploo、oracle。标准LOO为组内含截距CF分子与Chiang mHC3；p不计截距。CAR三层前10个系数分别为-0.4、0.5、0.2，层内计算并保留原层间项。样本Gram矩阵和回归交叉矩统一用n归一化，PLOO使用投影权重H；oracle使用已知总体协方差，权重为X_i' Sigma^{-1} X_j/n。CLOO/PLOO/oracle的方差按I、II分别计算，oracle的第一项采用(1+H_ij^oracle)H_ij^oracle。全部H1使用HR=exp(-0.5)，70%事件截止；非CAR试验1–10000，CAR试验20001–30000。完整设置见 `code/settings.csv`。

`raw/` 中 `n100_CRE_1to1_H0.rds` 对应Gaussian CRE、n=100、1:1、H0；BRE/2to1/H1类推。`n200/300/500/1000_CRE_1to1_*.rds` 为样本量实验，`gamma_n100_CRE_*.rds` 为Gamma实验，`sign_low_*.rds` 为CAR。

RDS采用xz无损压缩，可直接用 `readRDS()` 读取。`raw[试验,p,字段]` 保存各方法的u分子、v方差及部分标准化统计量；`parts[试验,p,分项]` 保存LOO分项。完整X和生存资料由 `experiment.R` 的种子规则重放。

`pooled_components[试验,p,方法,分项]` 保存CLOO/PLOO/oracle的I、II和between；总方差为三项之和。CRE/BRE的between为0，CAR保留层间项。对应均值汇总见 `pooled_variance_components.csv`。

`summary.csv` 一行对应group × p × method，`rate` 是拒绝率（0–1），`valid/invalid` 是实际有效/失效次数。BRE 2:1高维保留1536个未定义mHC3方差；拒绝率按各方法有效次数计算。每图CSV的 `value` 是画出的拒绝率。

## 使用

从本目录运行。R需要MASS和survival；Python需要numpy、pandas、matplotlib。

```powershell
# 绘制九张图
python code/plot.py

# 从保存raw核对全部700行汇总
Rscript --vanilla code/summarize.R

# 重放保存的记录
Rscript --vanilla code/run.R --mode=replay

# 运行某一组；已有同名raw时跳过，job见settings.csv
Rscript --vanilla code/run.R --mode=run --job=1 --B=10000

# 全部24组完成后，重新收集summary
Rscript --vanilla code/summarize.R --collect

# ACTG175 Z3应用：复算六方法及论文结果表
Rscript --vanilla application/code/run.R
```

实际应用的筛选条件、变量编码和输出字段见 [application/README.md](application/README.md)。数据已包含在仓库中，复算不需要连接外部数据服务。
