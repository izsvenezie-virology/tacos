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


def main(coverage_file: TextIO, output_file: str, min_coverage: int) -> None:
    """Plots the coverage data."""
    cov_df = pd.read_csv(
        coverage_file,
        sep="\t",
        header=None,
        names=["chrom", "pos", "cov"],
        keep_default_na=False,
    )

    fig = plt.figure(figsize=(16, 9))
    fig.subplots_adjust(left=0.05, right=0.98)
    fig = plot_coverage(fig, cov_df, min_coverage)

    fig.savefig(output_file)


def plot_coverage(fig: Figure, cov_df: DataFrame, min_coverage: int) -> Figure:
    """Plot the coverage split by chromosome with proportional subplot widths."""
    chroms = cov_df.chrom.unique().tolist()
    chrom_dfs = [cov_df.loc[cov_df.chrom == chrom] for chrom in chroms]

    width_ratios = [chrom_df.shape[0] for chrom_df in chrom_dfs]
    gs = fig.add_gridspec(
        2 if min_coverage > 0 else 1,
        len(chroms),
        width_ratios=width_ratios,
        hspace=0.25,
        wspace=0,
    )

    top_axes = plot_chroms_coverage(fig, gs, 0, chrom_dfs, 0)
    top_axes = format_top_axes(top_axes)
    if min_coverage > 0:
        bottom_axes = plot_chroms_coverage(fig, gs, 1, chrom_dfs, min_coverage)
        bottom_axes = format_bottom_axes(bottom_axes, min_coverage * 2)

    fig.supxlabel("Position", fontsize=16)

    return fig


def plot_chroms_coverage(
    fig: Figure,
    gs: GridSpec,
    row_idx: int,
    chrom_dfs: List[DataFrame],
    min_coverage: int,
) -> List[Axes]:
    """Plot chromosome coverage in a row of proportional subplots."""
    cmap = plt.get_cmap("Dark2")
    row_axes: List[Axes] = []

    for col_idx, chrom_df in enumerate(chrom_dfs):
        color = cmap(col_idx % cmap.N)

        ax = fig.add_subplot(
            gs[row_idx, col_idx], sharey=row_axes[0] if row_axes else None
        )

        ax.plot(chrom_df["pos"], chrom_df["cov"], linewidth=1, color=color)
        highlight_low_coverage_regions(ax, chrom_df, min_coverage)
        ax.set_title(chrom_df["chrom"].iloc[0], fontsize=8)
        ax = format_x_axis(ax, chrom_df.shape[0])
        ax = format_y_axis(ax)

        row_axes.append(ax)

    format_row_axes(row_axes)
    return row_axes


def format_y_axis(axes: Axes) -> Axes:
    axes.grid(which="major", alpha=0.8, color="#CCCCCC", linestyle="--")

    axes.spines["left"].set_visible(False)
    axes.spines["right"].set_visible(False)
    axes.tick_params(labelleft=False)
    axes.tick_params(axis="y", which="both", left=False, right=False)

    return axes


def highlight_low_coverage_regions(
    ax: Axes, chrom_df: DataFrame, min_coverage: int
) -> Axes:
    x = chrom_df["pos"].to_numpy()
    below_threshold = chrom_df["cov"].to_numpy() < min_coverage

    if below_threshold.any():
        ax.fill_between(
            x,
            0,
            1,
            where=below_threshold.tolist(),
            transform=ax.get_xaxis_transform(),
            color="#f6bcbc",
            alpha=0.35,
            interpolate=True,
            zorder=0,
        )
    return ax


def format_x_axis(ax: Axes, x_max: int) -> Axes:
    """Format X and Y axes parameters."""
    locator = mticker.MaxNLocator(nbins=4, min_n_ticks=3, integer=True)
    ticks = [int(round(tick)) for tick in locator.tick_values(0, x_max)]
    if 0 not in ticks:
        ticks.insert(0, 0)
    if x_max not in ticks:
        ticks.append(x_max)
    ticks = sorted(set(ticks))

    ax.set_xticks(ticks)
    tick_labels = [f"{tick:.0f}" for tick in ticks]
    tick_labels[0] = ""
    tick_labels[-2] = f"{x_max}-0"
    ax.set_xticklabels(tick_labels)
    ax.tick_params(axis="both", labelsize=8)
    ax.tick_params(axis="x", labelrotation=45)
    plt.setp(ax.get_xticklabels(), ha="right")
    ax.set_xlim(0, x_max)
    return ax


def format_row_axes(axes: List[Axes]) -> List[Axes]:
    first_ax = axes[0]
    first_ax.set_ylabel("Coverage", fontsize=16)
    first_ax.spines["left"].set_visible(True)
    first_ax.tick_params(axis="y", which="both", left=True, right=False, labelleft=True)
    tick_labels = first_ax.get_xticklabels()
    tick_labels[0].set_text("0")
    first_ax.set_xticklabels(tick_labels)

    last_ax = axes[-1]
    last_ax.spines["right"].set_visible(True)
    tick_labels = last_ax.get_xticklabels()
    tick_labels[-2].set_text(tick_labels[-2].get_text().replace("-0", ""))
    last_ax.set_xticklabels(tick_labels)
    return axes


def format_top_axes(axes: List[Axes]) -> List[Axes]:
    first_ax = axes[0]
    first_ax.set_yscale("log")
    return axes


def format_bottom_axes(axes: List[Axes], upper_lim: int) -> List[Axes]:
    last_ax = axes[-1]
    last_ax.set_ylim(0, upper_lim)
    return axes


# fmt: off
@click.command(help='Creates plots from a coverage file')
@click.version_option(__version__, '-v', '--version', message=f'%(prog)s, version %(version)s, by {__author__} ({__contact__})')
@click.help_option('-h', '--help')
@click.option('-m', '--min-coverage', default=0, show_default=True, help='Coverage threshold for marking low-coverage regions')
@click.argument('coverage_file', type=File('r'))
@click.argument('output_file', type=Path(file_okay=True))
def cli(*args, **kwargs):
    main(*args, **kwargs)


if __name__ == '__main__':
    cli()
