from typing import List, TextIO

import click
import matplotlib.ticker as mticker
import pandas as pd
from click.types import File, Path
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from pandas.core.frame import DataFrame

# Stops annoying warning messages form pandas
pd.options.mode.chained_assignment = None

__version__ = "2.0.0"
__author__ = "EdoardoGiussani"
__contact__ = "egiussani@izsvenezie.it"


def main(coverage_file: TextIO, output_file: str) -> None:
    """Plots the coverage data."""
    cov_df = pd.read_csv(
        coverage_file,
        sep="\t",
        header=None,
        names=["chrom", "pos", "cov"],
        keep_default_na=False,
    )

    fig = plt.figure(figsize=(16, 9), tight_layout=True)
    fig = plot_coverage(fig, cov_df)

    fig.savefig(output_file)


def plot_coverage(fig: Figure, cov_df: DataFrame) -> Figure:
    """Plot the coverage split by chromosome with proportional subplot widths."""
    chroms = cov_df.chrom.unique().tolist()
    chrom_dfs = [cov_df.loc[cov_df.chrom == chrom] for chrom in chroms]
    width_ratios = [chrom_df.shape[0] for chrom_df in chrom_dfs]

    gs = fig.add_gridspec(
        2,
        len(chroms),
        width_ratios=width_ratios,
        hspace=0.25,
        wspace=0,
    )

    top_axes = plot_chroms_coverage(fig, gs, 0, chroms, chrom_dfs)
    top_axes = format_subplot(top_axes, True)
    # top_axes[0].set_title("Sequence coverage", fontsize=20, pad=14)

    low_axes = plot_chroms_coverage(fig, gs, 1, chroms, chrom_dfs)
    low_axes = format_subplot(low_axes, False)
    for ax in low_axes:
        ax.set_ylim(top=100)

    fig.supxlabel("Position", fontsize=16)

    return fig


def plot_chroms_coverage(
    fig: Figure,
    gs: GridSpec,
    row_idx: int,
    chroms: List[str],
    chrom_dfs: List[DataFrame],
) -> List[Axes]:
    """Plot chromosome coverage in a row of proportional subplots."""
    row_axes: List[Axes] = []
    for col_idx, (chrom, chrom_df) in enumerate(zip(chroms, chrom_dfs)):
        ax = fig.add_subplot(
            gs[row_idx, col_idx], sharey=row_axes[0] if row_axes else None
        )
        ax.plot(chrom_df["pos"], chrom_df["cov"], linewidth=0.5)
        ax.set_title(chrom, fontsize=8)
        ax = format_axes(ax, chrom_df.shape[0])
        if row_axes:
            ax.tick_params(labelleft=False)
        row_axes.append(ax)

    return row_axes


def format_subplot(axes: List[Axes], log: bool) -> List[Axes]:
    """Apply common Y formatting to a row of axes."""
    last_idx = len(axes) - 1
    for idx, ax in enumerate(axes):
        if log:
            ax.set_yscale("log")
        ax.grid(which="major", alpha=0.8, color="#CCCCCC", linestyle="--")

        ax.spines["left"].set_visible(idx == 0)
        ax.spines["right"].set_visible(idx == last_idx)
        ax.tick_params(axis="y", which="both", left=idx == 0, right=False)

    axes[0].set_ylabel("Coverage", fontsize=16)
    return axes


def format_axes(ax: Axes, length: int) -> Axes:
    """Format X and Y axes parameters."""
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=4, min_n_ticks=3))
    ax.tick_params(axis="both", labelsize=8)
    ax.tick_params(axis="x", labelrotation=45)
    plt.setp(ax.get_xticklabels(), ha="right")
    ax.set_xlim(0, length)
    return ax


# fmt: off
@click.command(help='Creates plots from a coverage file')
@click.version_option(__version__, '-v', '--version', message=f'%(prog)s, version %(version)s, by {__author__} ({__contact__})')
@click.help_option('-h', '--help')
@click.argument('coverage_file', type=File('r'))
@click.argument('output_file', type=Path(file_okay=True))
def cli(*args, **kwargs):
    main(*args, **kwargs)


if __name__ == '__main__':
    cli()
