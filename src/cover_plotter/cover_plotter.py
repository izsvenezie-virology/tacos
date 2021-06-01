from matplotlib import pyplot as plt, colors
import click
import pandas as pd
import math
# Stops annoying warning messages form pandas
pd.options.mode.chained_assignment = None


@click.command(help='Creates plots from a coverage file')
@click.option('-o', '--output', type=click.Path(file_okay=True),
              help='name of output pdf file')
@click.option('-i', '--incremental', is_flag=True,
              help='creates a plot with the coverage values on x axis and the number of postionos on y axis')
@click.argument('coverage_file', type=click.types.File('r'))
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

    ax, ticks = plot_chroms_coverage(ax, cov_df)
    ax = format_subplot(ax, len(cov_df.chrom.unique()))
    ax = format_axes(ax, ticks)

    return fig


def plot_chroms_coverage(ax, cov_df):
    pos_offset = 0
    tick_pos = list()
    chroms = cov_df.chrom.unique()

    get_ticks_pos.distance = cov_df.shape[0]/40
    format_ticks.ticks = cov_df['pos'].tolist()

    for chrom in chroms:
        chrom_df = cov_df.loc[cov_df.chrom == chrom]
        chrom_df['plot_pos'] = chrom_df['pos'] + pos_offset
        chrom_len = chrom_df.shape[0]

        tick_pos += get_ticks_pos(chrom_len, pos_offset)

        ax.plot(chrom_df['plot_pos'], chrom_df['cov'],
                label=chrom, linewidth=.5)
        pos_offset += chrom_len

    tick_pos.append(cov_df.shape[0] - 1)
    return ax, tick_pos


def format_subplot(ax, legend_cols):
    if legend_cols > 10:
        legend_cols = 10
    ax.set_title('Sequence coverage', fontsize=20)
    ax.set_xlabel('Position', fontsize=16)
    ax.set_ylabel('Coverage', fontsize=16)
    ax.legend(loc='upper center', frameon=True, ncol=legend_cols)
    ax.grid(which='major', alpha=.8, color='#CCCCCC', linestyle='--')
    return ax


def format_axes(ax, tick_pos):
    ax.set_ylim(bottom=0)
    plt.locator_params(axis="y", nbins=20)

    ax.set_xticks(tick_pos)
    ax.set_xticklabels([], rotation=40, ha='right')
    ax.xaxis.set_major_formatter(format_ticks)
    ax.tick_params(axis='both', labelsize=8)
    return ax


def format_ticks(tick_val, _):
    int_val = int(tick_val)
    if int_val in range(len(format_ticks.ticks)):
        val = format_ticks.ticks[int_val]
        if val == 1 and int_val > 0:
            last_val = format_ticks.ticks[int_val - 1]
            val = f'{last_val}-{val}'
        return val
    return ''


def get_ticks_pos(length, offset):
    number = int(length/get_ticks_pos.distance)
    delta = math.ceil(length/number)

    pos = [offset]
    for tick in range(delta, length, delta):
        pos.append(tick + offset)
    return pos


def plot_incremental(cov_df, fig):
    chroms = cov_df.chrom.unique()
    color_pallette = list(colors.TABLEAU_COLORS.keys())
    gs = fig.add_gridspec(3, 3)
    x = 0
    y = 0

    for i, chrom in enumerate(chroms):
        incr_df = cov_df.loc[cov_df.chrom == chrom, 'cov'].value_counts(
        ).sort_index(ascending=True).to_frame()
        ax = fig.add_subplot(gs[y, x])
        ax.plot(incr_df.index.values.tolist(), incr_df['cov'].cumsum(
        ).sort_index(ascending=False), color_pallette[i])
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
