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

- 只读文件里的第一张逐级表（`stage target_Pa sig_plate_Pa ...`）。后面的增量表是 PFC 自己算的，这里不用，全部从原始数据重新计算，作为独立复核。
- 应变增量按 ΔH/H_prev 计算，和 `oedo_fn.p3dat`、规范里 Es 的定义一致。Es1-2 = 100 kPa / ((H100 − H200)/H100)，H100、H200 在加载段上插值得到。
- 孔隙比 e = A·H/Vs − 1，其中颗粒体积 Vs = W/(g·ρs)，W 取文件头部的 bed weight，ρs 默认 2600 kg/m³，可用 `--rho` 修改。a1-2 = (e100 − e200)/0.1 MPa。自检：(1+e1)/a1-2 应等于 Es1-2。
- 两端任一级出现 `hitmax=1` 的增量记为不可靠（`reliable=0`），M–σ 图上画成空心点。
- M_ur 取 1600 kPa 卸到 100 kPa 的割线模量。

### 已验证

- 用合成数据（kPa、Pa、只有 H 三种格式）检查过，结果和理论值一致。
- soft、med 的实测结果：复算出的 Es1-2 与 PFC 输出相差不到 0.5%（6.01 vs 6.00 MPa，21.1 vs 21.0 MPa）；按 (1+e1)/a1-2 算的也一致。

### 待核实

JGJ 79-2012 表 4.2.7 的垫层压缩模量参考值（砂约 20–30 MPa，碎石、卵石约 30–50 MPa）是凭记忆给的，网上没查到原表。写进论文前请对照规范原文。
