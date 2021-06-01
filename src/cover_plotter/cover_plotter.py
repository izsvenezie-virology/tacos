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
    cov_df = pd.read_csv(coverage_file, sep='\t', header=None, names=[
                         'chrom', 'pos', 'cov'], keep_default_na=False)

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
    tick_distance = cov_df.shape[0]/40
    tick_pos = list()
    tick_label = list()
    chroms = cov_df.chrom.unique()
    format_ticks.ticks = cov_df['pos'].tolist()

    for chrom in chroms:
        chrom_df = cov_df.loc[cov_df.chrom == chrom]
        chrom_df['plot_pos'] = chrom_df['pos'] + pos_offset
        chrom_len = chrom_df.shape[0]

        t_pos, t_lbl = get_ticks(chrom_len, tick_distance, pos_offset)
        tick_pos += t_pos
        tick_label += t_lbl

        ax.plot(chrom_df['plot_pos'], chrom_df['cov'],
                label=chrom, linewidth=.5)
        pos_offset += chrom_len

    tick_pos.append(cov_df.shape[0] - 1)
    tick_label.append(chrom_len)

    ax.set_title('Sequence coverage', fontsize=20)
    ax.set_xlabel('Position', fontsize=16)
    ax.set_ylabel('Coverage', fontsize=16)

    ax.legend(loc='upper center', frameon=True, ncol=len(chroms))
    ax.set_ylim(bottom=0)

    ax.yaxis.set_major_locator(ticker.LinearLocator(20))

    ax.xaxis.set_major_locator(ticker.LinearLocator(40))
    ax.xaxis.set_major_formatter(format_ticks)
    ax.tick_params(labelrotation=270)
    ax.set_xticks(tick_pos)

    ax.grid(which='major', alpha=.8, color='#CCCCCC', linestyle='--')

    return fig

def format_ticks(tick_val, tick_pos):
    if int(tick_val) in range(len(format_ticks.ticks)):
        return format_ticks.ticks[int(tick_val)]
    return ''

def get_ticks(length, distance, offset):
    number = int(length/distance)
    delta = math.ceil(length/number)

    pos = [1 + offset]
    labels = [1]
    if hasattr(get_ticks, 'last_max'):
        labels = [f'{get_ticks.last_max}-1']
    get_ticks.last_max = length

    for tick in range(1 + delta, length, delta):
        pos.append(tick + offset)
        labels.append(tick)

    return pos, labels


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
