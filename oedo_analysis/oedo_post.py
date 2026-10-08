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

# 列名别名（小写、去掉下划线和空格后匹配）。PFC 端写出的列名以 oedo_fn.p3dat 为准，
# 这里尽量放宽，匹配不到时会报错并列出实际表头。
ALIASES = {
    "stage": ["stage", "step", "level", "lvl", "i", "no", "idx"],
    "sv": ["sv", "sigv", "sigmav", "svtop", "sigvtop", "sigmavtop", "stop", "sigtop",
           "verticalstress", "sz", "sigz", "sigmaz", "ptop", "target"],
    "sv_bot": ["svbot", "sigvbot", "sigmavbot", "sbot", "sigbot", "svbase", "sigbase"],
    "H": ["h", "height", "hc", "thick", "thickness"],
    "eps": ["eps", "epsv", "ev", "strain", "epsz", "ez", "epsilon"],
    "sh": ["sh", "sigh", "sigmah", "sigx", "sigmax", "lateralstress", "shside", "sigside"],
}


def norm(s):
    return re.sub(r"[\s_\-\(\)\[\]]|kpa|pa|mpa", "", s.lower())


def is_num(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def read_table(path):
    """返回 (header 列表, 数值矩阵)。表头取最后一个非数值行。"""
    header, rows = None, []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith(";"):
                continue
            toks = [t for t in re.split(r"[,\s\t]+", s.lstrip("#").strip()) if t]
            if toks and all(is_num(t) for t in toks):
                rows.append([float(t) for t in toks])
            elif toks:
                header = toks
    if not rows:
        raise ValueError(f"{path}: 没有数值行")
    ncol = max(len(r) for r in rows)
    rows = [r for r in rows if len(r) == ncol]
    if header is None or len(header) != ncol:
        raise ValueError(f"{path}: 表头列数({0 if header is None else len(header)}) "
                         f"与数据列数({ncol})不一致，表头={header}")
    return header, np.array(rows)


def map_columns(header):
    keys = [norm(h) for h in header]
    col = {}
    for name, al in ALIASES.items():
        for i, k in enumerate(keys):
            if k in al and i not in col.values():
                col[name] = i
                break
    return col


def to_kpa(x):
    """σ 单位判断: 最大值 >1e4 视为 Pa（1600 kPa 级别对应 1.6e6 Pa）。"""
    return x / 1e3 if np.nanmax(np.abs(x)) > 1e4 else x


def interp_on_branch(sv, eps, target):
    """在单调加载段上按 σ 线性插值 ε。"""
    order = np.argsort(sv)
    s, e = sv[order], eps[order]
    if target < s[0] or target > s[-1]:
        return math.nan
    return float(np.interp(target, s, e))


def analyse(path):
    header, a = read_table(path)
    col = map_columns(header)
    if "sv" not in col:
        raise ValueError(f"{path}: 找不到竖向应力列，表头={header}")
    sv = to_kpa(a[:, col["sv"]])
    sv_bot = to_kpa(a[:, col["sv_bot"]]) if "sv_bot" in col else None
    sh = to_kpa(a[:, col["sh"]]) if "sh" in col else None
    H = a[:, col["H"]] if "H" in col else None
    if "eps" in col:
        eps = a[:, col["eps"]].astype(float)
        if np.nanmax(np.abs(eps)) > 1.0:  # 百分数
            eps = eps / 100.0
    elif H is not None:
        eps = (H[0] - H) / H[0]  # 以第一行（贴实级）高度为基准
    else:
        raise ValueError(f"{path}: 既无应变列也无高度列，表头={header}")

    # 加载段: σ 第一次达到最大值之前（含）；之后为卸载段
    imax = int(np.argmax(sv))
    load = slice(0, imax + 1)
    sv_l, eps_l = sv[load], eps[load]
    if np.any(np.diff(sv_l) <= 0):
        print(f"[警告] {os.path.basename(path)}: 加载段 σ 非严格递增，请检查伺服是否收敛",
              file=sys.stderr)

    stages = []
    n = len(sv)
    for i in range(n):
        d = {"i": i, "sv": sv[i], "eps": eps[i]}
        if H is not None:
            d["H"] = H[i]
        if sh is not None:
            d["sh"] = sh[i]
        if sv_bot is not None:
            d["sv_bot"] = sv_bot[i]
            d["top_bot_diff_pct"] = 100 * (sv[i] - sv_bot[i]) / sv[i] if sv[i] else math.nan
        if i > 0:
            ds, de = sv[i] - sv[i - 1], eps[i] - eps[i - 1]
            d["branch"] = "load" if i <= imax else "unload"
            d["M_MPa"] = ds / de / 1e3 if de != 0 else math.nan
            if sh is not None and ds != 0:
                k0 = (sh[i] - sh[i - 1]) / ds
                nu = k0 / (1 + k0)
                d["K0"], d["nu"] = k0, nu
                if math.isfinite(d["M_MPa"]) and nu < 0.5:
                    d["E_MPa"] = d["M_MPa"] * (1 + nu) * (1 - 2 * nu) / (1 - nu)
        else:
            d["branch"] = "seat"
        stages.append(d)

    e100 = interp_on_branch(sv_l, eps_l, 100.0)
    e200 = interp_on_branch(sv_l, eps_l, 200.0)
    es12 = (200 - 100) / (e200 - e100) / 1e3 if (e200 - e100) > 0 else math.nan

    # 100–200 kPa 段 K0（若有 σ_h）
    k0_12 = nu_12 = E12 = math.nan
    if sh is not None:
        h100 = float(np.interp(100.0, sv_l, sh[load]))
        h200 = float(np.interp(200.0, sv_l, sh[load]))
        k0_12 = (h200 - h100) / 100.0
        nu_12 = k0_12 / (1 + k0_12)
        if math.isfinite(es12) and nu_12 < 0.5:
            E12 = es12 * (1 + nu_12) * (1 - 2 * nu_12) / (1 - nu_12)

    # 卸载回弹模量: 峰值点到卸载段最后一点的割线
    mur = math.nan
    if imax < n - 1:
        ds = sv[imax] - sv[-1]
        de = eps[imax] - eps[-1]
        mur = ds / de / 1e3 if de > 0 else math.nan
    eps_res = eps[-1] if imax < n - 1 else math.nan

    tag = re.sub(r"^oedo_|\.txt$", "", os.path.basename(path))
    summary = {
        "tag": tag, "sv_max_kPa": float(sv[imax]), "eps_max": float(eps[imax]),
        "Es1_2_MPa": es12, "K0_1_2": k0_12, "nu_1_2": nu_12, "E_1_2_MPa": E12,
        "M_ur_MPa": mur, "eps_residual": eps_res,
    }
    return tag, header, col, stages, summary


def fmt(x, p=3):
    return "—" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{p}g}"


def write_outputs(results, outdir):
    os.makedirs(outdir, exist_ok=True)
    keys = ["tag", "sv_max_kPa", "eps_max", "Es1_2_MPa", "K0_1_2", "nu_1_2",
            "E_1_2_MPa", "M_ur_MPa", "eps_residual"]
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

    lines = ["| 垫层 | σmax/kPa | εmax | Es1-2/MPa | K0(100–200) | ν | E(100–200)/MPa | M_ur/MPa | 残余ε |",
             "|---|---|---|---|---|---|---|---|---|"]
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
            plt.rcParams["font.family"] = font_manager.FontProperties(fname=fp).get_name()
            break
    plt.rcParams["axes.unicode_minus"] = False

    fig, ax = plt.subplots(figsize=(5.5, 4))
    for tag, _, _, st, _ in results:
        pts = [(math.sqrt(d["sv"] * st[d["i"] - 1]["sv"]), d["M_MPa"]) for d in st
               if d.get("branch") == "load" and math.isfinite(d.get("M_MPa", math.nan))]
        if pts:
            x, y = zip(*pts)
            ax.loglog(x, y, "o-", label=tag)
    ax.set_xlabel("σv（级间几何平均）/ kPa")
    ax.set_ylabel("切线约束模量 M / MPa")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="center right", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "oedo_M_sigma.png"), dpi=200)

    fig, ax = plt.subplots(figsize=(5.5, 4))
    for tag, _, _, st, _ in results:
        ax.semilogx([d["sv"] for d in st], [100 * d["eps"] for d in st], "o-", label=tag)
    ax.invert_yaxis()
    ax.set_xlabel("σv / kPa")
    ax.set_ylabel("竖向应变 ε / %")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "oedo_e_logp.png"), dpi=200)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("resdir")
    ap.add_argument("-o", "--outdir", default=None)
    args = ap.parse_args()
    files = sorted(glob.glob(os.path.join(args.resdir, "oedo_*.txt")))
    files = [f for f in files if not os.path.basename(f).startswith("oedo_log")]
    order = {"soft": 0, "med": 1, "stiff": 2}
    files.sort(key=lambda f: next((v for k, v in order.items() if f"_{k}" in f), 9))
    if not files:
        sys.exit(f"{args.resdir} 下没有 oedo_*.txt")
    results = []
    for f in files:
        tag, header, col, st, s = analyse(f)
        print(f"[{tag}] 表头={header}  识别列={col}  共 {len(st)} 级")
        results.append((tag, header, col, st, s))
    outdir = args.outdir or os.path.join(args.resdir, "post")
    write_outputs(results, outdir)
    plot(results, outdir)
    print(f"输出 -> {outdir}")


if __name__ == "__main__":
    main()
