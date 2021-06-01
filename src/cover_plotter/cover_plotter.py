from click.types import File
from matplotlib import pyplot as plt, ticker, colors
import click
import pandas as pd
import math
pd.options.mode.chained_assignment = None # Stops annoying warning messages form pandas


@click.command(help='Creates plots from a coverage file')
@click.option('-o', '--output', type=click.Path(file_okay=True),
              help='the output pdf file path')
@click.option('-i', '--incremental', is_flag=True,
              help='creates a plot with the coverage values on x axis and the number of postionos on y axis')
@click.argument('coverage_file', type=File('r'))
def main(coverage_file, output, incremental):

    cov_df = pd.read_csv(coverage_file, sep='\t', 
                        header=None, names=['chrom', 'pos', 'cov'], 
                        keep_default_na=False)

    fig = plt.figure(figsize=(16, 9), tight_layout=True)
    if incremental:
        fig = plot_incremental(cov_df, fig)
    else:
        fig = plot_coverage(cov_df, fig)

    if output:
        fig.savefig(output)
    else:
        plt.show()


def plot_coverage(cov_df, fig):
    ax = fig.add_subplot(1, 1, 1)

    pos_offset = 0
    tick_pos = list()
    chroms = cov_df.chrom.unique()

    get_ticks.distance = cov_df.shape[0]/40
    format_ticks.ticks = cov_df['pos'].tolist()

    for chrom in chroms:
        chrom_df = cov_df.loc[cov_df.chrom == chrom]
        chrom_df['plot_pos'] = chrom_df['pos'] + pos_offset
        chrom_len = chrom_df.shape[0]

        tick_pos += get_ticks(chrom_len, pos_offset)

        ax.plot(chrom_df['plot_pos'], chrom_df['cov'],
                label=chrom, linewidth=.5)
        pos_offset += chrom_len

    tick_pos.append(cov_df.shape[0] - 1)

    ax.set_title('Sequence coverage', fontsize=20)
    ax.set_xlabel('Position', fontsize=16)
    ax.set_ylabel('Coverage', fontsize=16)

    ax.legend(loc='upper center', frameon=True, ncol=len(chroms))
    ax.set_ylim(bottom=0)

    ax.yaxis.set_major_locator(ticker.LinearLocator(20))

    ax.xaxis.set_major_formatter(format_ticks)
    ax.set_xticks(tick_pos)
    ax.tick_params(axis='x', labelrotation=90)

    ax.grid(which='major', alpha=.8, color='#CCCCCC', linestyle='--')

    return fig

def format_ticks(tick_val, tick_pos):
    int_val = int(tick_val)
    if int_val in range(len(format_ticks.ticks)):
        val = format_ticks.ticks[int_val]
        if val == 1 and int_val > 1:
            last_val = format_ticks.ticks[int_val - 1]
            val = f'{last_val}-{val}'
        return val
    return ''

def get_ticks(length, offset):
    number = int(length/get_ticks.distance)
    delta = math.ceil(length/number)

    pos = [0 + offset]
    for tick in range(1 + delta, length, delta):
        pos.append(tick + offset)
    return pos


def plot_incremental(cov_df, fig):
    chroms = cov_df.chrom.unique()
    color_pallette = list(colors.TABLEAU_COLORS.keys())
    gs = fig.add_gridspec(3,3)
    x = 0
    y = 0

    for i, chrom in enumerate(chroms):
        incr_df = cov_df.loc[cov_df.chrom == chrom, 'cov'].value_counts().sort_index(ascending=True).to_frame()
        ax = fig.add_subplot(gs[y,x])
        ax.plot(incr_df.index.values.tolist(), incr_df['cov'].cumsum().sort_index(ascending=False), color_pallette[i])
        ax.set_title(chrom)
        ax.set(xlabel='Coverage', ylabel='# positions with at least coverage')
        ax.grid(which='major', alpha=.8, color='#CCCCCC', linestyle='--')

        x += 1
        if x >= 3:
            x = 0
            y += 1
    
    fig.suptitle('Incremental coverage', fontsize=20)
    return fig


if __name__ == '__main__':
    main()
