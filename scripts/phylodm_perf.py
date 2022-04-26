###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import argparse
import os
import platform
import random
import re
import subprocess
import sys
from multiprocessing import Pool

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from dendropy.simulate import treesim
from tqdm import tqdm

__version__ = '1.1.0'


def print_help():
    lines = [f'phylodm_perf v{__version__}',
             'usage: [data_dir] [n_taxa] [n_trials] [cpus]']
    print('\n'.join(lines))


def tree_worker(args):
    path, n_taxa = args
    tree = treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5,
                                    num_extant_tips=n_taxa, rng=random.Random(42))
    with open(path, 'w') as f:
        f.write(tree.as_string(schema='newick')[5:])
    return


def queue_worker(args):
    out_path, tree_path, method, n_taxa, cur_trial = args
    exe = 'gtime' if platform.system() == 'Darwin' else '/usr/bin/time'
    script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'phylodm_perf_worker.py')
    args = [exe, '--verbose', 'python', script_path, tree_path, method]

    proc = subprocess.Popen(args, stderr=subprocess.PIPE, encoding='utf-8')
    stdout, stderr = proc.communicate()

    if proc.returncode == 0:
        with open(out_path, 'w') as fh:
            fh.write(stderr)
    else:
        raise Exception(f'Non-zero return code: {stderr}')
    return


def generate_test_data(data_dir, n_taxa, n_trials, cpus):
    # Generate a tree for each number of taxa.
    queue_trees = list()
    for cur_n_taxa in n_taxa:
        cur_tree_path = os.path.join(data_dir, f'{cur_n_taxa}.tree')
        if not os.path.isfile(cur_tree_path):
            queue_trees.append((cur_tree_path, cur_n_taxa))

    # Generate workers for each number of taxa and trials.
    queue_jobs = list()
    for method in ('dendropy', 'phylodm'):
        for cur_n_taxa in n_taxa:
            cur_tree_path = os.path.join(data_dir, f'{cur_n_taxa}.tree')
            for cur_trial in range(n_trials):
                cur_job_path = os.path.join(data_dir, f'{method}_{cur_n_taxa}_{cur_trial}.txt')
                if not os.path.isfile(cur_job_path):
                    queue_jobs.append((cur_job_path, cur_tree_path, method, cur_n_taxa, cur_trial))

    print('Generating trees.')
    with Pool(processes=cpus) as pool:
        list(tqdm(pool.imap_unordered(tree_worker, queue_trees), total=len(queue_trees)))

    print('Collecting data.')
    with Pool(processes=cpus) as pool:
        list(tqdm(pool.imap_unordered(queue_worker, queue_jobs), total=len(queue_jobs)))

    return


def parse_proc_output(path):
    with open(path) as fh:
        output = fh.read()

    wall_time = re.search(r'.*Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.+)', output).group(1)
    max_mem = int(re.search(r'.*Maximum resident set size \(kbytes\): (\d+)', output).group(1))
    exit_code = int(re.search(r'.*Exit status: (\d)', output).group(1))

    if exit_code != 0:
        print(output)
        raise Exception('Non zero exit code!!!')

    split_time = wall_time.split(':')
    if len(split_time) == 3:
        hr_sec = int(split_time[0]) * 3600
        mm_sec = int(split_time[1]) * 60
        sec = float(split_time[2])
        tot_sec = hr_sec + mm_sec + sec
    else:
        mm_sec = int(split_time[0]) * 60
        sec = float(split_time[1])
        tot_sec = mm_sec + sec

    return tot_sec, max_mem


def plot_test_data(data_dir, n_taxa, n_trials):
    # Extract the data from the text files.
    data = list()
    for method in ('dendropy', 'phylodm'):
        for cur_n_taxa in n_taxa:
            for cur_trial in range(n_trials):
                cur_job_path = os.path.join(data_dir, f'{method}_{cur_n_taxa}_{cur_trial}.txt')
                cur_sec, cur_mem_kb = parse_proc_output(cur_job_path)
                data.append((method, cur_n_taxa, cur_trial, cur_sec, cur_mem_kb))

    df = pd.DataFrame(data)
    df.columns = ['method', 'n_taxa', 'trial', 'tot_sec', 'max_mem']

    df['tot_sec'] = df['tot_sec'] / 60  # min
    df['max_mem'] = df['max_mem'] / 1024 / 1024  # gb

    df_phylo = df[df['method'] == 'phylodm']
    df_dendro = df[df['method'] == 'dendropy']

    # PLOT SETTINGS
    # sns.set_style("darkgrid")
    n_boot = 1e4

    # Plot each series overthe range of data points
    plt.figure(figsize=(7, 5))
    sns.set_style("whitegrid")

    # sns.set(style="ticks", palette="colorblind")
    sns.lineplot(x="n_taxa", y="tot_sec", data=df_phylo, markers=True,
                label='PhyloDM')
    ax = sns.lineplot(x="n_taxa", y="tot_sec", data=df_dendro,   markers=True,
                     label='DendroPy', color="r")
    ax.set(xlabel='Number of taxa in tree', ylabel='PDM construction time (min)',
           title='DendroPy vs. PhyloDM PDM Construction Time', yscale='log')

    ylabels = [f'{x:g}' for x in ax.get_yticks()]
    ax.set_yticklabels(ylabels)

    xlabels = [f'{int(x)}' for x in ax.get_xticks()]
    ax.set_xticklabels(xlabels)

    ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10.0, subs='all'))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    plt.grid(True, which="both", ls="--", c='gray', alpha=0.6)

    plt.legend()
    plt.savefig(os.path.join(data_dir, 'plot_cpu.png'))
    plt.show()

    # MEMORY
    # Fit the memory data
    plt.figure(figsize=(7, 5))

    sns.set_style("whitegrid")
    sns.lineplot(x="n_taxa", y="max_mem", data=df_phylo,   markers=True,
                label='PhyloDM')
    ax = sns.lineplot(x="n_taxa", y="max_mem", data=df_dendro,   markers=True,
                     label='DendroPy', color="r")
    ax.set(xlabel='Number of taxa in tree', ylabel='Maximum memory used (GB)',
           title='DendroPy vs. PhyloDM PDM Maximum Memory Usage', yscale='log')

    ylabels = [f'{x:g}' for x in ax.get_yticks()]
    ax.set_yticklabels(ylabels)

    xlabels = [f'{int(x)}' for x in ax.get_xticks()]
    ax.set_xticklabels(xlabels)

    ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10.0, subs='all'))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    plt.grid(True, which="both", ls="--", c='gray', alpha=0.6)

    plt.legend()
    plt.savefig(os.path.join(data_dir, 'plot_memory.png'))
    plt.show()
    return


def main(args=None):
    parser = argparse.ArgumentParser(description='Collect performance metrics.')
    parser.add_argument('data_dir', type=str, help='path to write data')
    parser.add_argument('n_taxa', type=str, help='number of taxa in each test: '
                                                 '[10,100,200]')
    parser.add_argument('n_trials', type=int, help='the number of trials: [10]')
    parser.add_argument('cpus', type=int, help='number of CPUs to use: [1]')

    # Verify that a subparser was selected
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f'phylodm_perf v{__version__}')
        sys.exit(0)
    else:
        print(f'phylodm_perf v{__version__}')
        args = parser.parse_args(args)

        # Create the output directory if it doesn't exist.
        if not os.path.isdir(args.data_dir):
            os.makedirs(args.data_dir)
        else:
            print('Using existing directory.')

        n_taxa = [int(x) for x in args.n_taxa.split(',')]

        generate_test_data(args.data_dir, n_taxa, args.n_trials, args.cpus)

        print('Plotting test data')
        plot_test_data(args.data_dir, n_taxa, args.n_trials)

    # Done - no errors.
    print('Done')
    sys.exit(0)


if __name__ == "__main__":
    main()
