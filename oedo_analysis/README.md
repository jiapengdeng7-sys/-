# 垫层侧限压缩试验后处理

在本机 PFC 5.0 GUI 中运行 `oedo_run.p3dat`，把结果交到云端，由 `oedo_post.py` 统一算出 Es1-2、M–σ 曲线和 K0/ν。

## 本机（PFC GUI 命令行）

```
call 'E:\PFC课题\PFC-3D\不同球形-刚性垫层-不同姿态\球形-刚性垫层-不同姿态\球形-长宽0.8 厚0.3\oedo_run.p3dat'
```

注意：console 版恢复 GUI 生成的存档时会在 `fish005_64.dll` 崩溃，所以必须用 GUI 版运行。

## 把结果交到云端

把下面这些文件放进本仓库的 `oedo_results/` 并推送（`.p3sav` 体积大，不用传）：

- `oedo_soft_w16_wf.txt`、`oedo_med_w16_wf.txt`、`oedo_stiff_w16_wf.txt`
- `oedo_log.txt`（里面有每级的收敛情况，以及 soft 垫层反算出的 E_c）
- 如果方便，再加上 `oedo_fn.p3dat`、`oedo_run.p3dat`，用来核对列名和单位

## 后处理

```
pip install numpy matplotlib
python3 oedo_analysis/oedo_post.py oedo_results
```

输出在 `oedo_results/post/`：

| 文件 | 内容 |
|---|---|
| `summary.md` / `summary.csv` | 每个垫层的 Es1-2、K0、ν、E（100–200 kPa）、回弹模量 M_ur、残余应变 |
| `stages_<tag>.csv` | 逐级的 σv、ε、σh、切线 M、增量 K0、ν、E，以及压板和底板的应力差 |
| `oedo_M_sigma.png` | 切线约束模量 M 随 σv 的变化（双对数坐标） |
| `oedo_e_logp.png` | ε–lg σv 压缩曲线 |

### 计算口径

- Es1-2 = 100 kPa / (ε200 − ε100)，取加载段数据；σ 不正好落在某一级上时线性插值。
- 应变基准：文件里有 ε 列就直接用（PFC 端定义）；只有高度 H 时，取第一行（5 kPa 贴实级）的高度为 H0。
- σ 单位：最大值超过 1e4 的按 Pa 处理，否则按 kPa。
- ν 由增量 K0 算出：ν = K0/(1+K0)；E = M(1+ν)(1−2ν)/(1−ν)。
- M_ur 取峰值 1600 kPa 到卸载终点 100 kPa 之间的割线模量。
- 列名不是严格匹配，有一组别名：先认 σv（压板），其次 σv 底板、H、ε、σh。识别结果会打印出来，如果识别错了，在 `ALIASES` 里补上对应的列名。

### 已验证

用已知答案的合成数据做了检查（M = 20/50/150 MPa，K0 = 0.4，分别用 kPa、Pa、只有 H 三种格式），Es1-2、ν = 0.286、E、M_ur 全部和理论值一致。

### 待核实

JGJ 79-2012 表 4.2.7 的垫层压缩模量参考值（砂约 20–30 MPa，碎石、卵石约 30–50 MPa）是凭记忆给的，网上没查到原表。写进论文前请对照规范原文。
