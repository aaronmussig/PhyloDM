import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

PATH_CSV = '/tmp/results.csv'


def fit_data(x, y):
    fit_param = np.polyfit(x, y, 2)
    model_pred = np.polyval(fit_param, x)
    abs_error = model_pred - y
    r_squared = 1.0 - (np.var(abs_error) / np.var(y))
    return fit_param, r_squared


def eqn_as_string(params):
    coeffs = list()
    for n, c in enumerate(params):
        if n > 0:
            if c >= 0:
                sign = '+'
            else:
                sign = ''
        else:
            sign = ''
        ex = len(params) - 1 - n
        if ex == 0:
            ex = ''
        elif ex == 1:
            ex = 'x'
        else:
            ex = f'x^{ex}'

        exp_notation = f'{c:.2e}'
        exp_coef, exp_pwr = exp_notation.split('e')
        exp_pwr = str(int(exp_pwr))
        if exp_pwr in {'-1', '1', '0'}:
            exp_string = f'{c:.4f}'
        else:
            exp_string = str(exp_coef) + 'e^{' + exp_pwr + '}'
        coeffs.append(f'{sign}{exp_string}{ex}')
    return ' '.join(coeffs)


def format_legend(label, eqn, r2):
    str_eqn = eqn_as_string(eqn)

    return f'{label} (${str_eqn}$, $R^2={r2:.2f}$)'


def main():
    df = pd.read_csv(PATH_CSV)
    df['tot_sec'] = df['tot_sec'] / 60
    df['max_mem'] = df['max_mem'] / 1024


    df_phylo = df[df['method'] == 'phylodm']
    df_dendro = df[df['method'] == 'dendropy']

    # Fit the time data
    phylo_time_fit, phylo_time_r2 = fit_data(x=df_phylo['n_taxa'], y=df_phylo['tot_sec'])
    dendro_time_fit, dendro_time_r2 = fit_data(x=df_dendro['n_taxa'], y=df_dendro['tot_sec'])

    # Plot each series over the range of data points
    plt.figure(figsize=(5, 3))

    sns.set(style="ticks", palette="colorblind")
    sns.regplot(x="n_taxa", y="tot_sec", data=df_phylo, x_estimator=np.mean, order=2,
                label='PhyloDM')
    ax = sns.regplot(x="n_taxa", y="tot_sec", data=df_dendro, x_estimator=np.mean, order=2,
                     label='DendroPy')
    ax.set(xlabel='Number of taxa in tree', ylabel='PDM construction time (min)',
           title='DendroPy vs. PhyloDM PDM Construction Time')
           #,yscale='log')

    xlabels = [f'{int(x)}' for x in ax.get_xticks()]
    ax.set_xticklabels(xlabels)

    plt.legend()
    plt.show()

    # Fit the memory data
    phylo_mem_fit, phylo_mem_r2 = fit_data(x=df_phylo['n_taxa'], y=df_phylo['max_mem'])
    dendro_mem_fit, dendro_mem_r2 = fit_data(x=df_dendro['n_taxa'], y=df_dendro['max_mem'])

    plt.figure(figsize=(5, 3))

    sns.set(style="ticks", palette="colorblind")
    sns.regplot(x="n_taxa", y="max_mem", data=df_phylo, x_estimator=np.mean, order=2,
                label='PhyloDM')
    ax = sns.regplot(x="n_taxa", y="max_mem", data=df_dendro, x_estimator=np.mean, order=2,
                     label='DendroPy')
    ax.set(xlabel='Number of taxa in tree', ylabel='Maximum memory used (GB)',
           title='DendroPy vs. PhyloDM PDM Maximum Memory Usage')

    xlabels = [f'{int(x)}' for x in ax.get_xticks()]
    ax.set_xticklabels(xlabels)

    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
