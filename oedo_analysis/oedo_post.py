#!/usr/bin/env python3
"""垫层侧限压缩试验（PFC oedo_run.p3dat 输出）后处理。

用法:
    python3 oedo_post.py <结果目录> [-o <输出目录>]

读取 <结果目录> 下的 oedo_<tag>.txt（tag 例如 soft_w16_wf / med_w16_wf / stiff_w16_wf），
不依赖 PFC 里已算好的 M/K0/ν 列，而是从原始的 σ_v、ε（或 H）、σ_h 重新计算:
  - 每级切线约束模量 M = Δσ_v/Δε
  - Es1-2: σ_v 由 100 kPa 到 200 kPa 的割线压缩模量（不在级上时线性插值）
  - 增量 K0 = Δσ_h/Δσ_v，ν = K0/(1+K0)，E = M(1+ν)(1-2ν)/(1-ν)
  - 卸载回弹模量 M_ur
输出: summary.csv、summary.md、per-stage CSV、图 oedo_M_sigma.png / oedo_e_logp.png。
"""
import argparse
import csv
import glob
import math
import os
import re
import sys

import numpy as np

# 列名别名（小写、去掉下划线、空格和单位后匹配）。oedo_fn.p3dat 实际写出的表头是
#   stage target_Pa sig_plate_Pa sig_floor_Pa H_m eps sig_h_Pa cycles aratio hitmax
ALIASES = {
    "stage": ["stage", "step", "level", "lvl", "i", "no", "idx"],
    "tgt": ["target", "tgt"],
    "sv": ["sigplate", "sv", "sigv", "sigmav", "svtop", "sigvtop", "sigmavtop", "sigtop",
           "verticalstress", "sz", "sigz", "sigmaz", "ptop"],
    "sv_bot": ["sigfloor", "svbot", "sigvbot", "sigmavbot", "sigbot", "svbase", "sigbase"],
    "H": ["hm", "h", "height", "hc", "thick", "thickness"],
    "eps": ["eps", "epsv", "ev", "strain", "epsz", "ez", "epsilon"],
    "sh": ["sigh", "sh", "sigmah", "sigx", "sigmax", "lateralstress", "sigside"],
    "cycles": ["cycles", "cyc"],
    "aratio": ["aratio"],
    "hitmax": ["hitmax"],
}


def norm(s):
    return re.sub(r"[\s_\-\(\)\[\]]|kpa$|pa$|mpa$|_m$", "", s.lower())


def is_num(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def read_table(path):
    """返回 (表头, 数值矩阵, 文件头部文字)。只取第一张表: 表头行之后连续的数值行。"""
    header, rows, pre = None, [], []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith(";"):
                continue
            toks = [t for t in re.split(r"[,\s\t]+", s.lstrip("#").strip()) if t]
            if toks and all(is_num(t) for t in toks):
                rows.append([float(t) for t in toks])
            elif rows:
                break  # 第一张表结束（后面是增量表等）
            else:
                if header is not None:
                    pre.append(" ".join(header))
                header = toks
    if not rows:
        raise ValueError(f"{path}: 没有数值行")
    ncol = max(len(r) for r in rows)
    rows = [r for r in rows if len(r) == ncol]
    if header is None or len(header) != ncol:
        raise ValueError(f"{path}: 表头列数({0 if header is None else len(header)}) "
                         f"与数据列数({ncol})不一致，表头={header}")
    return header, np.array(rows), "\n".join(pre)


def map_columns(header):
    keys = [norm(h) for h in header]
    col = {}
    for name, al in ALIASES.items():
        for a in al:  # 按别名优先级找
            hit = [i for i, k in enumerate(keys) if k == a and i not in col.values()]
            if hit:
                col[name] = hit[0]
                break
    return col


def to_kpa(x):
    """σ 单位判断: 最大值 >1e4 视为 Pa（1600 kPa 级别对应 1.6e6 Pa）。"""
    return x / 1e3 if np.nanmax(np.abs(x)) > 1e4 else x


def interp_on_branch(sv, y, target):
    """在单调加载段上按 σ 线性插值。"""
    order = np.argsort(sv)
    s, e = sv[order], y[order]
    if target < s[0] or target > s[-1]:
        return math.nan
    return float(np.interp(target, s, e))


def header_value(pre, pattern):
    m = re.search(pattern, pre)
    return float(m.group(1)) if m else math.nan


def analyse(path, rho_s=2600.0):
    header, a, pre = read_table(path)
    col = map_columns(header)
    if "sv" not in col:
        raise ValueError(f"{path}: 找不到竖向应力列，表头={header}")
    sv = to_kpa(a[:, col["sv"]])
    sv_bot = to_kpa(a[:, col["sv_bot"]]) if "sv_bot" in col else None
    sh = to_kpa(a[:, col["sh"]]) if "sh" in col else None
    H = a[:, col["H"]] if "H" in col else None
    hit = a[:, col["hitmax"]] if "hitmax" in col else np.zeros(len(sv))
    if H is not None:
        eps = (H[0] - H) / H[0]  # 与 oedo_fn 一致: 以第 1 级高度为基准（仅用于作图）
    elif "eps" in col:
        eps = a[:, col["eps"]].astype(float)
        if np.nanmax(np.abs(eps)) > 1.0:  # 百分数
            eps = eps / 100.0
    else:
        raise ValueError(f"{path}: 既无应变列也无高度列，表头={header}")

    def deps(i, j):
        """i -> j 的应变增量。有 H 时按 ΔH/H_i（规范 Es 定义，与 oedo_fn 相同）。"""
        if H is not None:
            return (H[i] - H[j]) / H[i]
        return (eps[j] - eps[i]) / (1 - eps[i])

    # 孔隙比: 由床重 W、颗粒密度和盒子面积得出颗粒体积
    W = header_value(pre, r"bed weight W\s*=\s*([-\d.eE+]+)")
    L = header_value(pre, r"box\s*([\d.]+)\s*x")
    Vs = W / 9.81 / rho_s if math.isfinite(W) else math.nan
    A = L * L if math.isfinite(L) else math.nan

    # 加载段: σ 第一次达到最大值之前（含）；之后为卸载段
    imax = int(np.argmax(sv))
    load = slice(0, imax + 1)
    sv_l = sv[load]
    if np.any(np.diff(sv_l) <= 0):
        print(f"[警告] {os.path.basename(path)}: 加载段 σ 非严格递增，请检查伺服是否收敛",
              file=sys.stderr)

    stages = []
    n = len(sv)
    for i in range(n):
        d = {"i": i + 1, "sv": sv[i], "eps": eps[i], "hitmax": int(hit[i])}
        if H is not None:
            d["H"] = H[i]
            if math.isfinite(Vs):
                d["e"] = A * H[i] / Vs - 1
        if sh is not None:
            d["sh"] = sh[i]
        if sv_bot is not None:
            d["sv_bot"] = sv_bot[i]
            d["bot_minus_top_kPa"] = sv_bot[i] - sv[i]
        if i > 0:
            ds, de = sv[i] - sv[i - 1], deps(i - 1, i)
            d["branch"] = "load" if i <= imax else "unload"
            d["M_MPa"] = ds / de / 1e3 if de != 0 else math.nan
            d["reliable"] = int(hit[i] == 0 and hit[i - 1] == 0)
            if sh is not None and ds != 0:
                k0 = (sh[i] - sh[i - 1]) / ds
                nu = k0 / (1 + k0)
                d["K0"], d["nu"] = k0, nu
                if math.isfinite(d["M_MPa"]) and nu < 0.5:
                    d["E_MPa"] = d["M_MPa"] * (1 + nu) * (1 - 2 * nu) / (1 - nu)
        else:
            d["branch"] = "seat"
        stages.append(d)

    # Es1-2: H(100), H(200) 在加载段上插值，Es = Δσ / (ΔH/H100)
    if H is not None:
        h1, h2 = (interp_on_branch(sv_l, H[load], t) for t in (100.0, 200.0))
        es12 = 100.0 / ((h1 - h2) / h1) / 1e3 if h1 > h2 else math.nan
    else:
        e1, e2 = (interp_on_branch(sv_l, eps[load], t) for t in (100.0, 200.0))
        es12 = 100.0 / ((e2 - e1) / (1 - e1)) / 1e3 if e2 > e1 else math.nan

    # 孔隙比口径的 a1-2 与 Es=(1+e1)/a1-2（应与上面的 Es1-2 相同，作为自检）
    e_100 = e_200 = a12 = math.nan
    if H is not None and math.isfinite(Vs):
        e_100, e_200 = A * h1 / Vs - 1, A * h2 / Vs - 1
        a12 = (e_100 - e_200) / 0.1  # MPa^-1
    # 100–200 kPa 段 K0（若有 σ_h）
    k0_12 = nu_12 = E12 = math.nan
    if sh is not None:
        h100 = interp_on_branch(sv_l, sh[load], 100.0)
        h200 = interp_on_branch(sv_l, sh[load], 200.0)
        k0_12 = (h200 - h100) / 100.0
        nu_12 = k0_12 / (1 + k0_12)
        if math.isfinite(es12) and nu_12 < 0.5:
            E12 = es12 * (1 + nu_12) * (1 - 2 * nu_12) / (1 - nu_12)

    # 卸载回弹模量: 峰值点到卸载段最后一点的割线；回弹比 = 卸载回弹量 / 峰值时总压缩量
    mur = math.nan
    if imax < n - 1:
        ds = sv[imax] - sv[-1]
        de = -deps(imax, n - 1)
        mur = ds / de / 1e3 if de > 0 else math.nan

    tag = re.sub(r"^oedo_|\.txt$", "", os.path.basename(path))
    summary = {
        "tag": tag, "n_hitmax": int(hit.sum()),
        "first_converged_kPa": float(sv[int(np.argmin(hit))]) if (hit == 0).any() else math.nan,
        "e_100": e_100, "e_200": e_200, "a1_2_perMPa": a12,
        "Es1_2_MPa": es12, "K0_1_2": k0_12, "nu_1_2": nu_12, "E_1_2_MPa": E12,
        "M_ur_MPa": mur,
    }
    return tag, header, col, stages, summary


def fmt(x, p=3):
    return "—" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{p}g}"


def write_outputs(results, outdir):
    os.makedirs(outdir, exist_ok=True)
    keys = ["tag", "n_hitmax", "first_converged_kPa", "e_100", "e_200", "a1_2_perMPa",
            "Es1_2_MPa", "K0_1_2", "nu_1_2", "E_1_2_MPa", "M_ur_MPa"]
    with open(os.path.join(outdir, "summary.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for _, _, _, _, s in results:
            w.writerow(s)
    for tag, _, _, stages, _ in results:
        sk = sorted({k for d in stages for k in d}, key=lambda k: (k != "i", k))
        with open(os.path.join(outdir, f"stages_{tag}.csv"), "w", newline="",
                  encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=sk)
            w.writeheader()
            w.writerows(stages)

    lines = ["| 垫层 | 未收敛级数 | 首个收敛级/kPa | e(100) | e(200) | a1-2/MPa⁻¹ | Es1-2/MPa "
             "| K0(100–200) | ν | E(100–200)/MPa | M_ur/MPa |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, _, _, _, s in results:
        lines.append("| " + " | ".join([s["tag"]] + [fmt(s[k]) for k in keys[1:]]) + " |")
    with open(os.path.join(outdir, "summary.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


def plot(results, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for fp in ["/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc",
               "C:/Windows/Fonts/simhei.ttf", "C:/Windows/Fonts/msyh.ttc"]:
        if os.path.exists(fp):
            font_manager.fontManager.addfont(fp)
            plt.rcParams["font.family"] = [font_manager.FontProperties(fname=fp).get_name(),
                                           "DejaVu Sans"]
            break
    plt.rcParams["axes.unicode_minus"] = False
    from matplotlib.ticker import FuncFormatter
    plain = FuncFormatter(lambda v, _: f"{v:g}")  # 避免 10^-1 的数学字体负号缺字

    colors = {"soft": "#1f77b4", "med": "#ff7f0e", "stiff": "#2ca02c"}

    def color(tag):
        return next((c for k, c in colors.items() if tag.startswith(k)), None)

    # M–σ: 级间几何平均应力；实心=两端都收敛，空心=含 hitmax 级（含压板空行程，不代表材料）
    fig, ax = plt.subplots(figsize=(5.5, 4))
    for tag, _, _, st, _ in results:
        c = color(tag)
        pts = [(math.sqrt(st[j]["sv"] * st[j - 1]["sv"]), st[j]["M_MPa"], st[j]["reliable"])
               for j in range(1, len(st))
               if st[j]["branch"] == "load" and math.isfinite(st[j]["M_MPa"])]
        if not pts:
            continue
        x, y, r = zip(*pts)
        ax.loglog(x, y, "-", color=c, alpha=0.5)
        ax.loglog([p[0] for p in pts if p[2]], [p[1] for p in pts if p[2]], "o",
                  color=c, label=tag)
        ax.loglog([p[0] for p in pts if not p[2]], [p[1] for p in pts if not p[2]], "o",
                  mfc="none", color=c)
    ax.plot([], [], "o", mfc="none", color="gray", label="含未收敛级（仅供参考）")
    ax.set_xlabel("σv（级间几何平均）/ kPa")
    ax.set_ylabel("切线约束模量 M / MPa")
    ax.xaxis.set_major_formatter(plain)
    ax.yaxis.set_major_formatter(plain)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "oedo_M_sigma.png"), dpi=200)

    # e–lg p 压缩曲线（无孔隙比时画 ε）
    fig, ax = plt.subplots(figsize=(5.5, 4))
    use_e = all("e" in st[0] for _, _, _, st, _ in results)
    for tag, _, _, st, _ in results:
        c = color(tag)
        y = [d["e"] if use_e else 100 * d["eps"] for d in st]
        x = [d["sv"] for d in st]
        k = next(i for i, d in enumerate(st) if d["branch"] == "unload") if any(
            d["branch"] == "unload" for d in st) else len(st)
        ax.semilogx(x[:k], y[:k], "o-", color=c, label=tag)
        ax.semilogx(x[k - 1:], y[k - 1:], "s--", color=c, mfc="none")
        for d, xi, yi in zip(st, x, y):
            if d["hitmax"]:
                ax.plot(xi, yi, "x", color="k", ms=8)
    ax.plot([], [], "x", color="k", label="hitmax（未收敛）")
    ax.plot([], [], "s--", color="gray", mfc="none", label="卸载")
    if not use_e:
        ax.invert_yaxis()
    ax.set_xlabel("σv / kPa")
    ax.set_ylabel("孔隙比 e" if use_e else "竖向应变 ε / %")
    ax.xaxis.set_major_formatter(plain)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "oedo_e_logp.png"), dpi=200)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("resdir")
    ap.add_argument("-o", "--outdir", default=None)
    ap.add_argument("--rho", type=float, default=2600.0, help="颗粒密度 kg/m3（cush_dens）")
    args = ap.parse_args()
    files = sorted(glob.glob(os.path.join(args.resdir, "oedo_*.txt")))
    files = [f for f in files if not os.path.basename(f).startswith("oedo_log")]
    order = {"soft": 0, "med": 1, "stiff": 2}
    files.sort(key=lambda f: next((v for k, v in order.items() if f"_{k}" in f), 9))
    if not files:
        sys.exit(f"{args.resdir} 下没有 oedo_*.txt")
    results = []
    for f in files:
        tag, header, col, st, s = analyse(f, args.rho)
        print(f"[{tag}] 表头={header}  识别列={col}  共 {len(st)} 级")
        results.append((tag, header, col, st, s))
    outdir = args.outdir or os.path.join(args.resdir, "post")
    write_outputs(results, outdir)
    plot(results, outdir)
    print(f"输出 -> {outdir}")


if __name__ == "__main__":
    main()
