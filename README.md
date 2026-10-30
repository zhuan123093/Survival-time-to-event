# 当前九图：代码与结果

按当前九张可见图整理。此仓库保存对应的代码与结果；论文和原项目文件保留在本地原位置。

```text
latest_9_figures_20261030/
├── README.md                 # 本文件：代码用途、数据说明、运行命令
├── INDEX.csv                 # 九图逐一对应图片、绘图CSV、原始结果
├── MANIFEST.csv              # 文件大小、SHA256及来源
├── code/
│   ├── settings.csv          # 24组设置：n、CRE/BRE/CAR、H0/H1、比例、p、试验偏移
│   ├── data.R                # 70%事件删失、风险集派生响应、标准log-rank
│   ├── methods.R             # 当前OLS/dof、CLOO/PLOO/oracle及共用矩阵运算
│   ├── chiang_loo.R          # 当前含截距CF分子、Chiang mHC3、饱和组及大样本计算
│   ├── experiment.R          # 当前资料生成、随机种子、七方法单次计算
│   ├── run.R                 # 按设置运行，或重放已保存记录核对代码
│   ├── summarize.R           # 从24份raw收集汇总，或核对保存的summary
│   └── plot.py               # 统一绘制当前九图，导出各图CSV
└── results/
    ├── raw/                  # 24份最终RDS，每组10000次，原文件逐字节复制
    ├── summary.csv           # 当前700个方法设置的汇总
    ├── variance_components.csv # 当前LOO方差分项
    ├── invalid_trials.csv    # 当前七方法保留的失效记录
    ├── figures/              # 九套正式PDF和PNG
    ├── plot_data/            # 每图一份CSV，共1134个绘图点
    ├── code_check.csv        # 精简代码的24组、74条记录、七方法回放核对
    ├── summary_check.csv     # raw重计与summary逐列核对
    └── plot_check.csv        # 九图数据、PNG像素与正式PDF核对
```

## 代码对应关系

`run.R` → `experiment.R`：按 `settings.csv` 生成当前资料，然后调用 `data.R` 计算风险集响应，调用 `methods.R` 和 `chiang_loo.R` 计算七种方法。全部依赖在本目录内。

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

## 当前设置和数据

七种方法为lr、ols、dof、loo、cloo、ploo、oracle。标准LOO为组内含截距CF分子与Chiang mHC3；p不计截距。CAR三层前10个系数分别为-0.4、0.5、0.2，层内计算并保留原层间项。PLOO沿用当前W=H。全部H1使用HR=exp(-0.5)，70%事件截止；非CAR试验1–10000，CAR试验20001–30000。完整设置见 `code/settings.csv`。

`raw/` 中 `n100_CRE_1to1_H0.rds` 对应Gaussian CRE、n=100、1:1、H0；BRE/2to1/H1类推。`n200/300/500/1000_CRE_1to1_*.rds` 为样本量实验，`gamma_n100_CRE_*.rds` 为Gamma实验，`sign_low_*.rds` 为当前CAR。

RDS的 `raw[试验,p,字段]` 保存各方法的u分子、v方差及部分标准化统计量；`parts[试验,p,分项]` 保存当前LOO分项。完整X和生存资料由 `experiment.R` 的种子规则重放。原始RDS的历史来源字段作为存档元数据保留，当前代码读取本目录的raw和设置运行。

`summary.csv` 一行对应group × p × method，`rate` 是拒绝率（0–1），`valid/invalid` 是实际有效/失效次数。BRE 2:1高维保留1536个未定义mHC3方差；拒绝率按各方法有效次数计算。每图CSV的 `value` 是画出的拒绝率。

## 使用

从本目录运行。R需要MASS和survival；Python需要numpy、pandas、matplotlib。

```powershell
# 当前九图
python code/plot.py

# 从保存raw核对全部700行汇总
Rscript --vanilla code/summarize.R

# 按24组保存记录核对精简代码
Rscript --vanilla code/run.R --mode=replay

# 运行某一组；已有同名raw时跳过，job见settings.csv
Rscript --vanilla code/run.R --mode=run --job=1 --B=10000

# 全部24组完成后，重新收集summary
Rscript --vanilla code/summarize.R --collect
```

本次只整理代码并核对保存结果。74条代表/饱和记录的七方法分子、方差与保存值一致，最大相对误差2.56e-17；700行汇总逐列核对通过；九张重画PNG与原版逐像素一致，正式PDF和raw均保留原字节。
