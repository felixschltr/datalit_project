import numpy as np
import pandas as pd
import cycler
import matplotlib as mpl
import matplotlib.pyplot as plt

from .summarize import get_diffs, get_means

# color settings
n = 20
#color = plt.cm.tab20(np.linspace(0,1,n))
mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', plt.cm.tab20.colors)

def plot_means_vs_diffs(
    rdd_data : pd.DataFrame,
    ymin_rdd : int,
    ymax_rdd : int,
    coi_rdd : str,
    nrg_data : pd.DataFrame,
    ymin_nrg : int,
    ymax_nrg : int,
    coi_nrg : str,
    **kwargs
):
    fig, ax = plt.subplots(figsize=kwargs.get('figsize', (10,8)))
    # get differences in share of renewables in energy data
    diff_df = get_diffs(nrg_data, ymin=ymin_nrg, ymax=ymax_nrg, coi=coi_nrg)
    # get means in share of renewables in RD&D data
    mean_df = get_means(rdd_data, ymin=ymin_rdd, ymax=ymax_rdd, coi=coi_rdd)
    # merge into singel data frame
    total_df = diff_df.merge(mean_df, how='inner', on='country')
    # compute overall correlationplt.cm.tab20(np.linspace(0,1,n))
    pearson_r = np.corrcoef(total_df['mean_budget'], total_df['share_difference'])[0][1]
    x_corr = min(total_df['mean_budget'])
    y_corr = max(total_df['share_difference'])
    
    for country in total_df['country'].unique():
        sub_data = total_df[total_df['country'] == country]
        x = sub_data['mean_budget']
        y = sub_data['share_difference']
        ax.scatter(x, y, label=country)
        # add country annotations
        ax.text(
            x + kwargs.get('text_offset_x', - (len(country) * 5 / 20)),
            y + kwargs.get('text_offset_y', 0.4),
            country,
            fontsize = kwargs.get('text_fontsize', 9)
            )
    # add correlation annotation
    ax.text(
        x_corr + 0.05,
        y_corr - 0.05,
        f'Pearson $r$ = {round(pearson_r, 3)}',
        fontsize = 12
        )
    ax.set_title(
        kwargs.get('title', 'RD&D spending vs. share of energy from renewable sources')
    )
    ax.set_ylabel(
        kwargs.get('x_label', f'Difference in share of energy from renewable sources between {ymin_nrg} and {ymax_nrg}\n(in %)')
    )
    ax.set_xlabel(
        kwargs.get('y_label', f'Mean share of total Energy Technology RD&D Budget invested in\nrenewable energy sources between {ymin_rdd} and {ymax_rdd}\n(in %)')
    )
    ax.grid(axis=kwargs.get('gridaxis', 'both'))
    ax.set_axisbelow(True)

def plot_rdd(rdd_data : pd.DataFrame, **kwargs):
    """
    Create plot of RD&D data for each country with year on the x-axis and
    share of budget invested in renewable energy sources on y-axis.

    Parameters
    ----------
    rdd_data: data on RD&D spending per country

    Keyword arguments:
    avg_window: {int} year window for computing rolling average
        defaults to 1, i.e. no averaging
    figsize: {tuple of int} defaults to (10,8)
    selected: {list/iterable of strings} selection of countries to plot
    plot_non_avg {boolean}: when True, plots true signal in background
    alpha {scalar} : transparency of ture signal (if plotted)
    gridaxis {str} : defaults to 'y'

    """
    fig, ax = plt.subplots(figsize=kwargs.get('figsize', (10,8)))
    for country in kwargs.get('selected', rdd_data['country'].unique()):
        sub_data = (
            rdd_data[(rdd_data['country'] == country)]
            .sort_values(by='year')
        )
        x = sub_data['year']
        y = sub_data['share_budget_renewables']
        avg_window = kwargs.get('avg_window', 1)
        y_avg = sub_data['share_budget_renewables'].rolling(window=avg_window).mean()
        ax.plot(x, y_avg, '-o', label=country)
        if kwargs.get('plot_non_avg', False):
            ax.plot(x, y, ls='--',alpha=kwargs.get('alpha', 0.6))
    ax.set_title('Share of energy technology RD&D budget invested in renewable energy sources\nby country and  year')
    ax.set_xlabel('Year')
    ylab = 'Share of RD&D budget invested in renewable energy sources (in %)'
    if avg_window != 1:
        ylab += f'\naveraged over windows of {avg_window} years'
    ax.set_ylabel(ylab)
    ax.legend(bbox_to_anchor=(1.0,1.0), loc='upper left')
    ax.grid(axis=kwargs.get('gridaxis', 'y'))
    ax.set_axisbelow(True)