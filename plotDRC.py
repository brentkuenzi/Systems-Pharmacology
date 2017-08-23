import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib
import pandas as pd
from typing import Callable, Tuple, Dict, Optional, Union

def plot_dose_response(df,
                        drug_column='drug name', response_column='response', dose_column='log(dose)', class_column='class', size=100,
                        n_cols=5, ylim=(0, None), horizontal_padding=0.1, vertical_padding=0.3,
                        style='whitegrid', font_scale=1.2,
                        class_to_color=None):
    """Plot a grid of nice-looking dose-response curves."""
    sns.set(font_scale=font_scale)
    sns.set_style(style)
    plot = sns.FacetGrid(df, col=drug_column, col_wrap=n_cols, hue=class_column, palette=class_to_color)
    plot.map(plt.plot, dose_column, response_column)
    size = size if isinstance(size, float) or isinstance(size, int) else df[size]
    plot.map(plt.scatter, dose_column, response_column, s=size)
    plot.set(xlim=(0, None), ylim=ylim)
    plot.set_titles(col_template="{col_name}")
    plot.fig.subplots_adjust(wspace=horizontal_padding, hspace=vertical_padding)
    return plot
