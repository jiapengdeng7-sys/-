#!/usr/bin/env python3
"""冲击算例后处理: python3 impact_post.py <结果目录> [标签片段, 默认 bedford_a0_t0_v30]

读 break_history / result_summary / energy_check，按 soft/med/stiff 汇总:
  时间步(由接触前下落速度反算)、破碎、首次停止时的侵入深度、底板峰值力、
  重叠量、能量分配、后期持续下沉速度。输出 impact_summary.md 和时程图。
"""
import glob, os, re, sys
import numpy as np

def table(path):
    L = [l.split() for l in open(path).read().splitlines() if l.strip()]
    return L[0], L[1:]

def main():
    rd = sys.argv[1]
    key = sys.argv[2] if len(sys.argv) > 2 else "bedford_a0_t0_v30"
    rows, hist = [], {}
    for tag in ["soft", "med", "stiff"]:
        bh = glob.glob(os.path.join(rd, f"break_history_cush{tag}_w16_wf_*{key}*.txt"))
        if not bh:
            continue
        h, b = table(bh[0]); c = {k: i for i, k in enumerate(h)}
        d = np.array([[float(x) for x in r] for r in b])
        _, es = table(bh[0].replace("break_history_", "energy_check_")); eh, _ = table(bh[0].replace("break_history_", "energy_check_"))
        e = dict(zip(eh, es[0]))
        v_in = float(e["v_in"]); m = float(e["m_cluster"]); Ein = float(e["E_in"])
        s = d[:, 0]; pen = d[:, c["penetration"]]; com = d[:, c["com_y"]]
        dt = (pen[2] - pen[1]) / (s[2] - s[1]) / v_in          # 接触前自由下落
        t = s * dt
        v = np.diff(com) / np.diff(t)
        i0 = int(np.argmax(pen > 0.04)); ia = i0 + int(np.argmax(v[i0:] >= 0)) + 1
        jf = int(np.argmax(d[:, c["fbase_y"]])); ko = int(np.argmax(d[:, c["ovl_max_r"]]))
        late = t >= 0.08
        sink = np.polyfit(t[late], com[late], 1)[0]
        ov = d[:, c["ovl_max_r"]]
        rows.append(dict(
            tag=tag, dt=dt, T=t[-1], n_bond0=int(d[0, c["n_intact"]]), crack_max=int(d[:, c["n_cracks"]].max()),
            crack_end=int(d[-1, c["n_cracks"]]), epb_max=d[:, c["epbstrain"]].max(),
            t_arrest=t[ia] * 1e3, pen_arrest=pen[ia], pen_end=pen[-1], clu_min_end=d[-1, c["clu_min_y"]],
            F_static=d[0, c["fbase_y"]], F_peak=d[jf, c["fbase_y"]], t_F=t[jf] * 1e3,
            ovl_max=ov[ko], t_ovl=t[ko] * 1e3, n10=int(d[:, c["n_ovl10"]].max()),
            ms_ovl50=(ov > 0.5).sum() * (s[1] - s[0]) * dt * 1e3,
            slip=float(e["eslip"]) / Ein, dash=float(e["edashpot"]) / Ein, ek_end=float(e["ekinetic"]) / Ein,
            body=float(e["ebody"]), resid=float(e["resid_over_Ein"]), sink=sink,
            v_clu_end=np.sqrt(2 * float(e["E_kin_cluster"]) / m)))
        hist[tag] = (t, d, c)
    hdr = ["垫层", "破碎键数(最大/末)", "pb应变能峰值/J", "首次停止 t/ms", "首次停止侵入/m", "末态侵入/m",
           "底板峰值力/kN (静载)", "峰值时刻/ms", "最大重叠/r_min", ">0.5r 持续/ms", "n(>10%)峰值",
           "摩擦耗能", "阻尼耗能", "末动能", "能量残差", "后期下沉/(m/s)"]
    out = ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    for r in rows:
        out.append("| " + " | ".join([
            r["tag"], f"{r['crack_max']}/{r['crack_end']} (共{r['n_bond0']})", f"{r['epb_max']:.3g}",
            f"{r['t_arrest']:.1f}", f"{r['pen_arrest']:.3f}", f"{r['pen_end']:.3f}",
            f"{r['F_peak']/1e3:.1f} ({r['F_static']/1e3:.1f})", f"{r['t_F']:.1f}",
            f"{r['ovl_max']:.2f}", f"{r['ms_ovl50']:.1f}", f"{r['n10']}",
            f"{r['slip']:.1%}", f"{r['dash']:.1%}", f"{r['ek_end']:.1%}", f"{r['resid']:.1e}",
            f"{-r['sink']:.2f}"]) + " |")
    txt = "\n".join(out)
    print(f"dt = {rows[0]['dt']:.3e} s, 总时长 {rows[0]['T']:.3f} s\n" + txt)
    open(os.path.join(rd, "impact_summary.md"), "w", encoding="utf-8").write(txt + "\n")

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib import font_manager
    fp = "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc"
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
        plt.rcParams["font.family"] = [font_manager.FontProperties(fname=fp).get_name(), "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    col = {"soft": "#1f77b4", "med": "#ff7f0e", "stiff": "#2ca02c"}
    fig, ax = plt.subplots(2, 2, figsize=(10, 7))
    for tag, (t, d, c) in hist.items():
        ms = t * 1e3
        ax[0, 0].plot(ms, d[:, c["fbase_y"]] / 1e3, color=col[tag], label=tag)
        ax[0, 1].plot(ms, d[:, c["penetration"]], color=col[tag], label=tag)
        ax[1, 0].plot(ms, d[:, c["ovl_max_r"]], color=col[tag], label=tag)
        ax[1, 1].plot(ms, d[:, c["ekinetic"]], color=col[tag], label=tag)
    ax[0, 0].set(xlim=(0, 40), xlabel="t / ms", ylabel="底板力 / kN", title="底板竖向力（前 40 ms）")
    ax[0, 1].set(xlabel="t / ms", ylabel="侵入深度 / m（相对最高颗粒）", title="前缘侵入深度（全程）")
    ax[1, 0].axhline(0.1, color="gray", ls="--", lw=0.8); ax[1, 0].axhline(1.0, color="red", ls="--", lw=0.8)
    ax[1, 0].set(xlim=(0, 40), xlabel="t / ms", ylabel="最大重叠 / r_min", title="垫层相关接触最大重叠（前 40 ms）")
    ax[1, 1].set(yscale="log", xlabel="t / ms", ylabel="全部颗粒动能 / J", title="系统动能（全程）")
    for a in ax.flat:
        a.grid(alpha=0.3); a.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(rd, "impact_history.png"), dpi=160)

if __name__ == "__main__":
    main()
