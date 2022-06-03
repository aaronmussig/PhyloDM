import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer

plt.rcParams['font.sans-serif'] = "Arial"


def plot_joint_axis(path_results, path_out):
    df = pd.read_csv(path_results, sep="\t")

    df = df[(df['n_taxa'] >= 1000)]
    df['tot_min'] = df['tot_sec'] / 60
    df['max_mem_gb'] = df['max_mem'] / 1024 / 1024

    df_dendro = df[df['source'] == 'dendropy']
    df_phylodm = df[df['source'] == 'phylodm']
    df_loading = df[df['source'] == 'loading']

    phylodm_cpu = df_phylodm[['n_taxa', 'tot_min']].groupby("n_taxa").median()
    dendropy_cpu = df_dendro[['n_taxa', 'tot_min']].groupby("n_taxa").median()
    loading_cpu = df_loading[['n_taxa', 'tot_min']].groupby("n_taxa").median()

    phylodm_mem = df_phylodm[['n_taxa', 'max_mem_gb']].groupby("n_taxa").median()
    dendropy_mem = df_dendro[['n_taxa', 'max_mem_gb']].groupby("n_taxa").median()
    loading_mem = df_loading[['n_taxa', 'max_mem_gb']].groupby("n_taxa").median()

    t = list(phylodm_cpu.index)
    phylodm_cpu_vals = [float(max(0, x)) for x in phylodm_cpu['tot_min'].values - loading_cpu['tot_min'].values]
    dendropy_cpu_vals = [float(max(0, x)) for x in dendropy_cpu['tot_min'].values - loading_cpu['tot_min'].values]

    phylodm_mem_vals = [float(max(0, x)) for x in phylodm_mem['max_mem_gb'].values - loading_mem['max_mem_gb'].values]
    dendropy_mem_vals = [float(max(0, x)) for x in dendropy_mem['max_mem_gb'].values - loading_mem['max_mem_gb'].values]

    fig, ax1 = plt.subplots(figsize=(8, 6))

    color_dendro = '#3727B3'
    color_phylo = '#179A39'

    ax1.set_xlabel('n_taxa')
    ax1.set_ylabel('Time (minutes) / Memory (GB)')
    p1 = ax1.plot(t, phylodm_cpu_vals, color=color_phylo, label='PhyloDM (minutes)')
    p2 = ax1.plot(t, dendropy_cpu_vals, color=color_dendro, label='DendroPy (minutes)')
    ax1.grid()
    ax1.set_yscale('log')

    p3 = ax1.plot(t, phylodm_mem_vals, color=color_phylo, label='PhyloDM (GB)', linestyle='dashed')
    p4 = ax1.plot(t, dendropy_mem_vals, color=color_dendro, label='DendroPy (GB)', linestyle='dashed')

    ax1.grid(True, which="both", ls="--", c='gray', alpha=0.6)

    ylabels = [f'{x:g}' for x in ax1.get_yticks()]
    ax1.set_yticklabels(ylabels)

    font_size = 12

    # Set maximum values
    ax1_phylo_max_y = max(phylodm_cpu_vals)
    ax1_phylo_max_x = max(t)
    ax1.text(ax1_phylo_max_x, ax1_phylo_max_y, f'{ax1_phylo_max_y * 60:.1f} seconds',
             fontsize=font_size, color=color_phylo, va='bottom', ha='right')

    ax1_dendro_max_y = max(dendropy_cpu_vals)
    ax1_dendro_max_x = max(t)
    ax1.text(ax1_dendro_max_x, ax1_dendro_max_y - 50, f'{ax1_dendro_max_y / 60:.1f} hours',
             fontsize=font_size, color=color_dendro, va='top', ha='right')

    # Set maximum values
    ax2_phylo_max_y = max(phylodm_mem_vals)
    ax2_phylo_max_x = max(t)
    ax1.text(ax2_phylo_max_x, ax2_phylo_max_y, f'{ax2_phylo_max_y:.1f} GB',
             fontsize=font_size, color=color_phylo, va='bottom', ha='right')

    ax2_dendro_max_y = max(dendropy_mem_vals)
    ax2_dendro_max_x = max(t)
    ax1.text(ax2_dendro_max_x, ax2_dendro_max_y, f'{ax2_dendro_max_y:.1f} GB',
             fontsize=font_size, color=color_dendro, va='bottom', ha='right')

    # fig.tight_layout()

    plots = p1 + p3 + p2 + p4
    labs = [p.get_label() for p in plots]
    ax1.legend(plots, labs, loc=0)

    plt.title('PhyloDM vs DendroPy distance matrix construction resource usage')

    # Polyfit
    print_eqn(t, phylodm_cpu_vals)
    print_eqn(t, dendropy_cpu_vals)
    print_eqn(t, dendropy_mem_vals)
    print_eqn(t, phylodm_mem_vals)

    plt.savefig(path_out)


def print_eqn(x, y):
    cpu_polyfit = np.polyfit(x, y, 2)
    cpu_f = np.poly1d(cpu_polyfit)
    print(f'y = {cpu_f.c[0]} x^2 + {cpu_f.c[1]} x + {cpu_f.c[2]}')


def main(path_results: str, path_out: str):
    # plot_subplot(path_results, path_out)
    plot_joint_axis(path_results, path_out)


if __name__ == "__main__":
    typer.run(main)
