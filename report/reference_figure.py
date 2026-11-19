"""Render the report figures from the distributed chr22 scores and aggregate results."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator

from pldsc import meta_ld_scores, meta_weights, single_reference_bias


HERE = Path(__file__).resolve().parent
INK, TEAL, RUST, GREY, GRID = "#1E3641", "#007668", "#A6472B", "#596C77", "#DCE4E6"
GROUPS = ("AFR", "EAS", "SAS", "AMR")


def frame(ax, title, grid_axis="y"):
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold", color=INK, pad=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#9AA9AF")
    ax.tick_params(length=3, color="#8A9CA4")
    ax.grid(axis=grid_axis, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def save(fig, name):
    fig.savefig(HERE / name, dpi=300, facecolor="white")
    plt.close(fig)


def reference_figure():
    with np.load(HERE.parent / "examples/mini/chr22_ldscores.npz") as data:
        names = list(data["ancestries"])
        eur = names.index("EUR")
        scores = data["scores"]
        available = data["available"] & (data["maf"] >= float(data["maf_gwas"]))
        keep = data["maf"][eur] >= float(data["maf_regression"])
    reference = scores[eur, eur, keep]
    rho = np.full((len(names), len(names)), 0.8)
    np.fill_diagonal(rho, 1.0)
    shares = (0.10, 0.25, 0.50)
    residuals = np.empty((len(GROUPS), len(shares)))
    fig, (scatter, matrix) = plt.subplots(1, 2, figsize=(7.25, 3.1), gridspec_kw={"width_ratios": [1.04, 1]})
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.18, top=0.86, wspace=0.40)
    print("Other group   Share   OLS kappa   Residual score variance (%)")
    for i, name in enumerate(GROUPS):
        for j, share in enumerate(shares):
            fractions = np.zeros(len(names))
            fractions[eur], fractions[names.index(name)] = 1 - share, share
            omega, n = meta_weights(fractions, available)
            matched = meta_ld_scores(scores, omega, rho=rho)[keep]
            n = n[keep]
            x, y = n * reference / n.mean(), n * matched / n.mean()
            slope = single_reference_bias(matched, reference, sample_sizes=n).kappa
            residuals[i, j] = 100 * (1 - np.corrcoef(x, y)[0, 1]**2)
            print(f"{name:>11}   {share:.2f}       {slope:.3f}       {residuals[i, j]:.2f}")
            if name == "AFR" and share == 0.25:
                scatter.scatter(x, y, s=3.5, alpha=0.16, color=TEAL, linewidths=0)
                grid = np.array([0, x.max()])
                scatter.plot(grid, grid, color=GREY, lw=1.25, ls="--", label="Equal scores")
                scatter.plot(grid, y.mean() + slope * (grid - x.mean()), color=RUST, lw=1.5,
                             label=f"Fitted slope {slope:.3f}")
    frame(scatter, "A   75% EUR and 25% AFR")
    scatter.set(xlim=(0, 64), ylim=(0, 64), xticks=[0, 20, 40, 60], yticks=[0, 20, 40, 60],
                xlabel="EUR score, scaled by sample size", ylabel="Matched score, same scaling")
    scatter.legend(loc="upper left", frameon=False, fontsize=8, handlelength=2)
    frame(matrix, "B   Variation beyond EUR")
    matrix.grid(False)
    cmap = LinearSegmentedColormap.from_list("residual", ["#F0F6F3", "#93C9BD", "#146C65"])
    matrix.imshow(residuals, cmap=cmap, vmin=0, vmax=9, aspect="auto", interpolation="nearest")
    matrix.set(xticks=range(3), xticklabels=["10%", "25%", "50%"], yticks=range(4),
               yticklabels=GROUPS, xlabel="Other group's share of sample size")
    matrix.spines[:].set_visible(False)
    matrix.tick_params(length=0)
    for i in range(len(GROUPS)):
        for j in range(len(shares)):
            value = residuals[i, j]
            matrix.text(j, i, f"{value:.2f}%", ha="center", va="center",
                        color="white" if value > 7.5 else INK, fontsize=10,
                        fontweight="bold" if (i, j) == (0, 1) else "normal")
    matrix.add_patch(Rectangle((0.52, -0.48), 0.96, 0.96, fill=False, color=RUST, lw=1.6))
    save(fig, "reference_composition.png")


def asthma_figure():
    data = json.loads((HERE / "asthma_score_shape.json").read_text())
    stats, rows = data["summary"], data["chromosomes"]
    critical = 2.0796138447  # 95% t with 21 df
    fig, (scale, risk) = plt.subplots(1, 2, figsize=(7.25, 3.05), gridspec_kw={"width_ratios": [0.88, 1.36]})
    fig.subplots_adjust(left=0.10, right=0.985, bottom=0.21, top=0.84, wspace=0.40)
    frame(scale, "A   Asthma heritability", "x")
    for y, (key, color) in enumerate([("matched_h2", TEAL), ("eur_h2", GREY)]):
        value, se = stats[key]["estimate"], stats[key]["se"]
        scale.errorbar(value, y, xerr=critical * se, fmt="o", color=color, lw=1.7, capsize=3, ms=6)
        scale.annotate(f"{value:.3f}", (value, y), xytext=(0, 8), textcoords="offset points",
                       ha="center", va="bottom", fontsize=9, fontweight="bold", color=color)
    scale.set(yticks=[0, 1], yticklabels=["Matched", "EUR"], ylim=(-0.52, 1.53),
              xlabel="Balanced observed-scale SNP $h^2$")
    scale.xaxis.set_major_locator(MaxNLocator(nbins=4))
    scale.tick_params(axis="y", length=0)
    scale.spines["left"].set_visible(False)
    frame(risk, "B   Held-out prediction error")
    chrom = [r["chr"] for r in rows]
    values = np.array([100 * (r["eur_sse"] - r["matched_sse"]) / r["eur_sse"] for r in rows])
    risk.axhline(0, color=GREY, lw=0.95)
    risk.vlines(chrom, 0, values, colors=np.where(values >= 0, TEAL, RUST), lw=1)
    risk.scatter(chrom, values, c=np.where(values >= 0, TEAL, RUST), s=12, zorder=3)
    gain = stats["gain_percent"]
    risk.axvline(23.3, color=GRID, lw=1)
    risk.errorbar(25, gain["estimate"], yerr=critical * gain["se"], color=INK,
                  fmt="D", ms=5, capsize=3.3, lw=1.3)
    risk.set(xticks=[1, 5, 10, 15, 20, 25], xticklabels=["1", "5", "10", "15", "20", "All"],
             xlim=(0, 26.6), xlabel="Chromosome", ylabel="Reduction in EUR error (%)")
    risk.yaxis.set_major_locator(MaxNLocator(nbins=5))
    save(fig, "asthma_prediction.png")


def wilson_interval(hits, total):
    p, z = hits / total, 1.9599639845
    divisor = 1 + z*z/total
    center = (p + z*z/(2*total)) / divisor
    radius = z * np.sqrt(p*(1-p)/total + z*z/(4*total**2)) / divisor
    return 100 * np.array([center - radius, center + radius])


def simulation_figure():
    data = json.loads((HERE / "score_shape_simulation.json").read_text())
    n = data["replicates"]
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 3.0))
    fig.subplots_adjust(left=0.13, right=0.98, bottom=0.29, top=0.83, wspace=0.35)
    styles = [(None, "Population LD", INK, "D", 0.18), (400, "400 people", RUST, "s", 0),
              (6000, "6,000 people", TEAL, "o", -0.18)]
    for ax, model, title, target, limit in zip(axes, ("EUR", "matched"),
            ("A   False positives", "B   Power"), (5, 80), ((0, 10), (0, 103))):
        frame(ax, title, "x")
        ax.axvline(target, ls="--", color=GREY, lw=1)
        for panel_n, label, color, marker, offset in styles:
            cells = sorted((r for r in data["cells"] if r["panel_n"] == panel_n and r["generating_model"] == model),
                           key=lambda r: r["other_share"])
            hits = np.array([r["rejected_zero"] for r in cells])
            p = 100 * hits / n
            lo, hi = wilson_interval(hits, n)
            ax.errorbar(p, np.array([2, 1, 0]) + offset, xerr=[p - lo, hi - p], color=color,
                        fmt=marker, capsize=2, ms=4.3, lw=1.1, label=label)
        ax.set(ylim=(-0.45, 2.48), xlim=limit, yticks=[2, 1, 0],
               yticklabels=["90 : 10", "75 : 25", "50 : 50"] if model == "EUR" else [],
               xlabel="E-only null rejected (%)" if model == "EUR" else "Matched pattern detected (%)")
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
    axes[0].set_ylabel("Simulated E : F mixture")
    axes[0].set_xticks([0, 2.5, 5, 7.5, 10], ["0", "2.5", "5", "7.5", "10"])
    axes[1].set_xticks([0, 20, 40, 60, 80, 100])
    fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", ncol=3,
               frameon=False, fontsize=8.5, title="LD reference per group", title_fontsize=8.5,
               bbox_to_anchor=(0.51, -0.002))
    save(fig, "score_shape_power.png")


def main():
    plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans", "text.color": INK,
                         "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK})
    reference_figure()
    asthma_figure()
    simulation_figure()


if __name__ == "__main__":
    main()
