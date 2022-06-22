from bokeh.palettes import Set3

def _get_palette(cmap=Set3[12], n=12, start=0, end=1):
    import matplotlib
    import numpy as np

    linspace = np.linspace(start, end, n)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("customcmap", cmap)
    palette = cmap(linspace)
    hex_palette = sorted([matplotlib.colors.rgb2hex(c) for c in palette])
    return hex_palette

