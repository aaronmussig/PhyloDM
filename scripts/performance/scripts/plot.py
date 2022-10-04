import matplotlib.patheffects as PathEffects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from scipy.optimize import curve_fit

plt.rcParams['font.sans-serif'] = "Arial"


def plot_joint_axis(path_results, path_out):
    df = pd.read_csv(path_results, sep="\t")

    df = df[(df['n_taxa'] >= 1000)]
    df['tot_min'] = df['tot_sec'] / 60
    df['max_mem_gb'] = df['max_mem'] / 1024 / 1024

    df_dendro = df[df['source'] == 'dendropy']
    df_phylodm = df[df['source'] == 'phylodm']
    df_loading = df[df['source'] == 'loading']

    loading_cpu = df_loading[['n_taxa', 'tot_min']].groupby("n_taxa").median()
    loading_mem = df_loading[['n_taxa', 'max_mem_gb']].groupby("n_taxa").median()

    d_loading_cpu_mins = loading_cpu['tot_min'].to_dict()
    d_loading_mem_gb = loading_mem['max_mem_gb'].to_dict()

    t = list(map(int, df_phylodm['n_taxa'].values))
    phylodm_cpu_vals = [max(row['tot_min'] - d_loading_cpu_mins[row['n_taxa']], 0)
                        for idx, row in df_phylodm[['n_taxa', 'tot_min']].iterrows()]
    dendropy_cpu_vals = [max(row['tot_min'] - d_loading_cpu_mins[row['n_taxa']], 0)
                         for idx, row in df_dendro[['n_taxa', 'tot_min']].iterrows()]

    phylodm_mem_vals = [max(row['max_mem_gb'] - d_loading_mem_gb[row['n_taxa']], 0)
                        for idx, row in df_phylodm[['n_taxa', 'max_mem_gb']].iterrows()]
    dendropy_mem_vals = [max(row['max_mem_gb'] - d_loading_mem_gb[row['n_taxa']], 0)
                         for idx, row in df_dendro[['n_taxa', 'max_mem_gb']].iterrows()]

    fig, ax1 = plt.subplots(figsize=(8, 6))

    color_dendro = '#3727B3'
    color_phylo = '#179A39'
    markersize = 4

    # Polyfit
    popt_phylo_cpu, pcov = curve_fit(func, t, phylodm_cpu_vals)
    print(f'PhyloDM (cpu): x^2 * {popt_phylo_cpu[0]}')

    popt_phylo_mem, pcov = curve_fit(func, t, phylodm_mem_vals)
    print(f'PhyloDM (mem): x^2 * {popt_phylo_mem[0]}')

    popt_dendro_cpu, pcov = curve_fit(func, t, dendropy_cpu_vals)
    print(f'DendroPy (cpu): x^2 * {popt_dendro_cpu[0]}')

    popt_dendro_mem, pcov = curve_fit(func, t, dendropy_mem_vals)
    print(f'DendroPy (mem): x^2 * {popt_dendro_mem[0]}')

    alpha = 1

    xp = np.linspace(1000, 30000, 1000)
    p1 = ax1.plot(xp, func(xp, *popt_phylo_cpu), '-', color=color_phylo, label=f'PhyloDM (minutes) = x\u00b2 {popt_phylo_cpu[0]:.2e}')
    p2 = ax1.plot(xp, func(xp, *popt_phylo_mem), '--', color=color_phylo, label=f'PhyloDM (GB) = x\u00b2 {popt_phylo_mem[0]:.2e}')
    p3 = ax1.plot(xp, func(xp, *popt_dendro_cpu), '-', color=color_dendro, label=f'DendroPy (minutes) = x\u00b2 {popt_dendro_cpu[0]:.2e}')
    p4 = ax1.plot(xp, func(xp, *popt_dendro_mem), '--', color=color_dendro, label=f'DendroPy (GB) = x\u00b2 {popt_dendro_mem[0]:.2e}')

    ax1.set_xlabel('n_taxa')
    ax1.set_ylabel('Time (minutes) / Memory (GB)')
    ax1.plot(t, phylodm_cpu_vals, 'o', markersize=markersize, color=color_phylo, label='PhyloDM (minutes)', alpha=alpha)
    ax1.plot(t, dendropy_cpu_vals, 'o', markersize=markersize, color=color_dendro, label='DendroPy (minutes)', alpha=alpha)
    ax1.grid()
    ax1.set_yscale('log')

    ax1.plot(t, phylodm_mem_vals, 'x', markersize=markersize, color=color_phylo, label='PhyloDM (GB)', alpha=alpha)
    ax1.plot(t, dendropy_mem_vals, 'x', markersize=markersize, color=color_dendro, label='DendroPy (GB)', alpha=alpha)

    ax1.grid(True, which="both", ls="--", c='gray', alpha=0.6)

    ylabels = [f'{x:g}' for x in ax1.get_yticks()]
    ax1.set_yticklabels(ylabels)

    font_size = 12

    # Set maximum values
    ax1_phylo_max_x = max(t)
    ax1_phylo_max_y = func(ax1_phylo_max_x, popt_phylo_cpu)[0]
    txt = ax1.text(ax1_phylo_max_x, ax1_phylo_max_y, f'{ax1_phylo_max_y * 60:.1f} seconds',
                   fontsize=font_size, color=color_phylo, va='bottom', ha='right')
    txt.set_path_effects([PathEffects.withStroke(linewidth=4, foreground='w')])

    ax1_dendro_max_x = max(t)
    ax1_dendro_max_y = func(ax1_dendro_max_x, popt_dendro_cpu)[0]
    txt = ax1.text(ax1_dendro_max_x, ax1_dendro_max_y - 50, f'{ax1_dendro_max_y / 60:.1f} hours',
                   fontsize=font_size, color=color_dendro, va='top', ha='right')
    txt.set_path_effects([PathEffects.withStroke(linewidth=4, foreground='w')])

    # Set maximum values
    ax2_phylo_max_x = max(t)
    ax2_phylo_max_y = func(ax2_phylo_max_x, popt_phylo_mem)[0]
    txt = ax1.text(ax2_phylo_max_x, ax2_phylo_max_y, f'{ax2_phylo_max_y:.1f} GB',
                   fontsize=font_size, color=color_phylo, va='bottom', ha='right')
    txt.set_path_effects([PathEffects.withStroke(linewidth=4, foreground='w')])

    ax2_dendro_max_x = max(t)
    ax2_dendro_max_y = func(ax2_dendro_max_x, popt_dendro_mem)[0]
    txt = ax1.text(ax2_dendro_max_x, ax2_dendro_max_y, f'{ax2_dendro_max_y:.1f} GB',
                   fontsize=font_size, color=color_dendro, va='bottom', ha='right')
    txt.set_path_effects([PathEffects.withStroke(linewidth=4, foreground='w')])

    # fig.tight_layout()

    plots = p1 + p2 + p3 + p4
    labs = [p.get_label() for p in plots]
    ax1.legend(plots, labs, loc=0)

    plt.title('PhyloDM vs DendroPy distance matrix construction resource usage')
    plt.rcParams['svg.fonttype'] = 'none'

    # plt.show()
    plt.savefig(path_out)


def func(x, a):
    return x * x * a


def print_eqn(x, y):
    cpu_polyfit = np.polyfit(np.log(x), y, 1)
    # cpu_f = np.poly1d(cpu_polyfit)
    print(f'y = {cpu_polyfit[0]} log(x) + {cpu_polyfit[1]}')
    return cpu_polyfit[0], cpu_polyfit[1]


def main(path_results: str, path_out: str):
    # plot_subplot(path_results, path_out)
    plot_joint_axis(path_results, path_out)


if __name__ == "__main__":
    typer.run(main)
