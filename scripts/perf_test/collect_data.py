import os
import re
import subprocess
from multiprocessing import Pool

import pandas as pd
from tqdm import tqdm


def run_proc(n_taxa, method):
    script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'run_test.py')
    args = ['/usr/bin/time', '--verbose', '/srv/home/uqamussi/.conda/envs/py3/bin/python', script_path, str(n_taxa),
            method]
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            encoding='utf-8')
    return proc.communicate()


def parse_proc_output(output):
    wall_time = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.+)', output).group(1)
    max_mem = int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', output).group(1))
    exit_code = int(re.search(r'Exit status: (\d)', output).group(1))

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


def run_proc_worker(job):
    n_taxa, method, trial = job
    stdout, stderr = run_proc(n_taxa, method)
    tot_sec, max_mem = parse_proc_output(stderr)
    return {'n_taxa': n_taxa, 'method': method, 'tot_sec': tot_sec,
            'max_mem': max_mem / 1024, 'trial': trial}


def main():
    trials = range(10)
    points = [50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 4000, 6000, 8000, 10000]

    # Generate a queue of jobs
    queue = list()
    for n_taxa in points:
        for trial in trials:
            for method in ['phylodm', 'dendropy']:
                queue.append((n_taxa, method, trial))

    # Collect results
    pool = Pool(processes=40)
    results = list()
    for result in tqdm(pool.imap_unordered(run_proc_worker, queue), total=len(queue)):
        results.append(result)
    df = pd.DataFrame(results)
    df.to_csv('results.csv')


if __name__ == '__main__':
    main()
