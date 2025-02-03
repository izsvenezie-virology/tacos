import math
from typing import Any, List, Tuple

import click
import pandas as pd
from click.types import File, Path
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from pandas.core.frame import DataFrame

# Stops annoying warning messages form pandas
pd.options.mode.chained_assignment = None

__version__ = '2.0.0'
__author__ = 'EdoardoGiussani'
__contact__ = 'egiussani@izsvenezie.it'

@click.command(help='Creates plots from a coverage file')
@click.version_option(__version__, '-v', '--version', message=f'%(prog)s, version %(version)s, by {__author__} ({__contact__})')
@click.help_option('-h', '--help')
@click.option('-o', '--output', type=Path(file_okay=True), help='Name of output pdf file.')
@click.argument('coverage_file', type=File('r'))
def main(coverage_file: File, output: Path) -> None:
    '''Plots the coverage data.'''
    cov_df = pd.read_csv(coverage_file, sep='\t',
                         header=None, names=['chrom', 'pos', 'cov'],
                         keep_default_na=False)

    fig = plt.figure(figsize=(16, 9), tight_layout=True)
    fig = plot_coverage(cov_df, fig)

    if output:
        fig.savefig(output)
    else:
        plt.show()


def plot_coverage(cov_df: DataFrame, fig: Figure) -> Figure:
    '''Plot the coverage in standard way'''
    ax = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)

    ax, ticks = plot_chroms_coverage(ax, cov_df)
    ax = format_subplot(ax, len(cov_df.chrom.unique()))
    ax = format_axes(ax, ticks)

    ax2, ticks2 = plot_chroms_coverage(ax2, cov_df)
    ax2 = format_subplot(ax2, len(cov_df.chrom.unique()))
    ax2 = format_axes(ax2, ticks)
    ax2.set_ylim(top=100)
    ax2.legend().remove()

    return fig


def plot_chroms_coverage(ax: plt.axes, cov_df: DataFrame) -> Tuple[plt.axes, List[int]]:
    '''Plots the coverage for each chromosome and set x ticks'''
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


def format_subplot(ax:plt.axes, legend_cols: int) -> plt.axes:
    '''Format the axes setting title, labels, legend and grid'''
    if legend_cols > 10:
        legend_cols = 10
    ax.set_title('Sequence coverage', fontsize=20)
    ax.set_xlabel('Position', fontsize=16)
    ax.set_ylabel('Coverage', fontsize=16)
    ax.set_yscale('log')
    ax.legend(loc='upper center', frameon=True, ncol=legend_cols)
    ax.grid(which='major', alpha=.8, color='#CCCCCC', linestyle='--')
    return ax


def format_axes(ax: plt.axes, tick_pos: List[int]) -> plt.axes:
    '''Format X and Y axes parameters'''
    ax.set_xticks(tick_pos)
    ax.set_xticklabels([], rotation=40, ha='right')
    ax.xaxis.set_major_formatter(format_ticks)
    ax.tick_params(axis='both', labelsize=8)
    return ax


def format_ticks(tick_val: str, _: Any) -> str:
    '''Formatter for X ticks'''
    int_val = int(tick_val)
    if int_val in range(len(format_ticks.ticks)):
        val = format_ticks.ticks[int_val]
        if val == 1 and int_val > 0:
            last_val = format_ticks.ticks[int_val - 1]
            val = f'{last_val}-{val}'
        return val
    return ''


def get_ticks_pos(length: int, offset: int) -> List[int]:
    '''Collect all positions for X ticks'''
    number = int(length/get_ticks_pos.distance)
    if number == 0:
        number = 1
    delta = math.ceil(length/number)

    pos = [offset]
    for tick in range(delta, length, delta):
        pos.append(tick + offset)
    return pos


if __name__ == '__main__':
    main()
