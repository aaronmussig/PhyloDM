import os
import re

import pandas as pd
import typer


def parse_trial_output(path):
    with open(path) as fh:
        output = fh.read()

    wall_time = re.search(r'.*Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.+)', output).group(1)
    max_mem = int(re.search(r'.*Maximum resident set size \(kbytes\): (\d+)', output).group(1))
    exit_code = int(re.search(r'.*Exit status: (\d)', output).group(1))

    if exit_code != 0:
        print(output)
        raise Exception(f'Non zero exit code: {path}')

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


def get_trials_from_dir(dir_path: str, source: str):
    rows = list()
    for file in os.listdir(dir_path):
        if not file.endswith('.txt'):
            continue

        n_taxa, trial = file.split('_')
        n_taxa = int(n_taxa)
        trial = int(trial.split('.')[0])

        path = os.path.join(dir_path, file)
        tot_sec, max_mem = parse_trial_output(path)

        rows.append({
            'source': source,
            'n_taxa': n_taxa,
            'trial': trial,
            'tot_sec': tot_sec,
            'max_mem': max_mem
        })
    return pd.DataFrame(rows)


def main(dendropy_dir: str, phylodm_dir: str, loading_dir: str, output_path: str):
    df_dendropy = get_trials_from_dir(dendropy_dir, 'dendropy')
    df_phylodm = get_trials_from_dir(phylodm_dir, 'phylodm')
    df_loading = get_trials_from_dir(loading_dir, 'loading')

    df = pd.concat([df_phylodm, df_dendropy, df_loading], ignore_index=True)
    df.sort_values(by=['source', 'n_taxa', 'trial'], inplace=True)
    df.to_csv(output_path, index=False, sep='\t')


if __name__ == "__main__":
    typer.run(main)
