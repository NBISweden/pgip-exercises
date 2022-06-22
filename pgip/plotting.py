# NB: https://speciationgenomics.github.io/pca/ calculates explained
# variance from sums of plink eigenvals, although not all components
# have been calculated
import re
from bokeh.io import output_file
from bokeh.layouts import gridplot
from bokeh.plotting import figure

from bokeh.models import ColumnDataSource
import itertools
import math

def bokeh_plot_pca_coords(df, explained, *, pc1=1, pc2=2, **kw):
    """Plot pca coordinates with bokeh"""
    source = ColumnDataSource(df)
    x = explained[pc1 - 1]
    y = explained[pc2 - 1]
    xlab = f"PC{pc1} ({x:.1f}%)"
    ylab = f"PC{pc2} ({y:.1f}%)"

    metadata_columns = [x for x in df.columns if not re.match("^(PC[0-9]+|color)$", x)]
    tooltips = [(x, f"@{x}") for x in metadata_columns]
    p = figure(
        x_axis_label=xlab,
        y_axis_label=ylab,
        tooltips=tooltips,
        title=f"PC{pc1} vs PC{pc2}",
        **kw,
    )
    p.circle(
        x=f"PC{pc1}",
        y=f"PC{pc2}",
        source=source,
        color="color",
        size=15,
        alpha=0.8,
        line_color="black",
        legend_group="population",
    )
    p.add_layout(p.legend[0], "right")
    return p


def bokeh_plot_pca(df, eigenvals, ncomp=6, filename=None, **kw):
    """Make PCA plot with bokeh"""
    pairs = list(itertools.combinations(range(ncomp), 2))
    n = len(pairs)
    ncols = kw.pop("ncols", math.floor(math.sqrt(n)))
    plots = []
    for (i, j) in pairs:
        p = bokeh_plot_pca_coords(df, eigenvals, pc1=i + 1, pc2=j + 1, **kw)
        plots.append(p)
    gp = gridplot(plots, ncols=ncols)
    if filename is not None:
        output_file(filename)
        show(gp)
    else:
        return gp

