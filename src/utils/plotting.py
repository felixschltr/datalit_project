
from turtle import color
from typing import List
import numpy as np
import pandas as pd
import cycler
import matplotlib as mpl
import matplotlib.pyplot as plt
from textwrap import wrap
from tueplots import bundles
from .summarize import get_diffs, get_means

# tueplots settings
# mpl.rcParams.update(bundles.neurips2021(usetex=False))
# color settings
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
    if kwargs.get('correl', True):
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

def plot_lines(
    rdd_data : pd.DataFrame,
    nrg_data : pd.DataFrame,
    t_window = 10,
    y_shift = 10,
    avg_window = 1,
    **kwargs
    ):
    fig, ax1 = plt.subplots(figsize=kwargs.get('figsize', (10,8)))
    ax2 = ax1.twiny() # create twin y axis
    # ax1 for nrg_data
    ymax1 = kwargs.get('ymax1',max(nrg_data['year'].unique()))
    ymin1 = ymax1 - t_window
    assert ymin1 >= min(nrg_data['year'].unique()), \
        f"bounds of available data exceeded!"
    nrg_data = nrg_data[nrg_data['year'].isin(range(ymin1, ymax1+1))]
    # ax2 for rdd_data
    ymax2 = ymax1 - y_shift
    ymin2 = ymax2 - t_window
    assert ymin2 >= min(rdd_data['year'].unique()), \
        f"bounds of available data exceeded!"
    rdd_data = rdd_data[rdd_data['year'].isin(range(ymin2, ymax2+1))]
    ###
    for country in kwargs.get('selected', rdd_data['country'].unique()):
        # plot nrg_data
        label1 = 'share of primary energy from renewable sources'
        sub_data = (
            nrg_data[(nrg_data['country'] == country)]
            .sort_values(by='year')
        )
        x1 = sub_data['year']
        y1 = sub_data['prim_nrg_share_renewables']
        ax1.plot(x1, y1, '-s', c='C0', label=label1)
        # ax1 formatting
        x_ticks = [int(y) for y in x1]
        ax1.set_xticks(x_ticks)
        ax1.set_xticklabels(x_ticks, rotation=45)
        # plot rdd_data
        sub_data = (
            rdd_data[(rdd_data['country'] == country)]
            .sort_values(by='year')
        )
        x2 = sub_data['year']
        y2 = sub_data['share_budget_renewables']
        #avg_window = kwargs.get('avg_window', 1)
        y2_avg = y2.rolling(window=avg_window).mean()
        # plot non-averaged data
        label2 = 'share of RD&D spending on renewable energy technology'
        if avg_window != 1:
            label2 += ' (averaged)'
        if kwargs.get('plot_non_avg', False):
            ax2.plot(x2, y2, ls='--',alpha=kwargs.get('alpha', 0.6), c='C2', label='share RD&D non-average')
        # plot averaged data
        ax2.plot(x2, y2_avg, '-o', c='C2', label=label2)
        # ax2 formatting
        x_ticks = [int(year) for year in x2]
        ax2.set_xticks(x_ticks)
        ax2.set_xticklabels(x_ticks, rotation=45)
        # add correlation annotation
        if kwargs.get('correl', False):
            #get correlation
            pearson_r = np.corrcoef(y1,y2)[0][1]
            x_corr = min(x1)
            y_corr = max(max(y1), max(y2_avg))
            ax1.text(
                x_corr + 0.05,
                y_corr - 0.05,
                f'Pearson $r$ = {round(pearson_r, 3)}',
                fontsize = 12
                )
    # global formatting
    default_title = (
        'Share of RD&D budged invested in renewable energy sources\nvs.'+
        '\nShare of renewables in primary energy\ntime shifted'+
        f'by {y_shift} years'
    )
    ax1.set_xlabel(f'Time for {label1}')
    ax2.set_xlabel(f'Time for {label2}')
    ax1.set_ylabel('Share in %') # global y-axis label
    ax1.set_title(kwargs.get('title', default_title))
    lin1, lab1 = ax1.get_legend_handles_labels()
    lin2, lab2, = ax2.get_legend_handles_labels()
    labels = ['\n'.join(wrap(l, 10)) for l in lab1+lab2]
    ax1.legend(
        lin1 + lin2, labels, bbox_to_anchor=(1.0,1.0),
        loc='upper left', labelspacing=1)
    ax1.grid(axis=kwargs.get('gridaxsis', 'y'))
    ax1.set_axisbelow(True)

def get_mean_std(data : pd.DataFrame, selected : List[str], coi : str) -> pd.DataFrame:
    """
    Returns per-year means and (corrected) standard deviations of values
    in column of interest over selected countries

    Parameters
    ----------
    data : data frame containing data sorted by year and country
    selected : list of countries to compute means and std.s for
    coi : name of column of interest, i.e. column to compute statistics for

    Returns
    -------
    merged data frame containing mean and (corrected) st.dev. per year
    """
    data = data[data['country'].isin(selected)]
    mean_df = (
        data.groupby('year').mean()
        .reset_index()
        [['year', coi]]
    )
    std_df = (
        data.groupby('year').std(ddof=1)
        .reset_index()
        [['year', coi]]
    )
    # merge into single data frame
    merged_df = mean_df.merge(std_df, how='inner', on='year')
    # rename columns
    merged_df.rename(columns={f'{coi}_x':'mean', f'{coi}_y':'std'}, inplace=True)
    return merged_df

def plot_lines_global(
    rdd_data,
    nrg_data,
    t_window : int = 10,
    y_shift : int = 10,
    **kwargs
    ):
    # formatting of error bars
    elinewidt = kwargs.get('elinewidt', 1)
    capsize = kwargs.get('capsize', 2)
    # select countries
    selected = kwargs.get('selected', rdd_data['country'].unique())
    ## get statistics for rdd_data
    rdd_coi = kwargs.get('rdd_coi', 'share_budget_renewables')
    rdd_data = get_mean_std(rdd_data, selected, rdd_coi)
    ## get statistics for nrg_data
    nrg_coi = kwargs.get('nrg_coi','prim_nrg_share_renewables')
    nrg_data = get_mean_std(nrg_data, selected, nrg_coi)
    ## acutal plotting
    fig, ax1 = plt.subplots(figsize=kwargs.get('figsize', (10,8)))
    ax2 = ax1.twiny() # create twin y axis
    # ax1 for nrg_data
    label1 = 'share of primary energy from renewable sources'
    ymax1 = kwargs.get('ymax1',max(nrg_data['year'].unique()))
    ymin1 = ymax1 - t_window
    assert ymin1 >= min(nrg_data['year'].unique()), \
        f"bounds of available data exceeded!"
    nrg_data = nrg_data[nrg_data['year'].isin(range(ymin1, ymax1+1))]
    x1 = nrg_data['year']
    y1 = nrg_data['mean']
    y1err = nrg_data['std']
    ax1.errorbar(
        x1, y1, y1err, fmt='-s', elinewidth=elinewidt,
        capsize=capsize, color='C0', label=label1)
    # ax1 formatting
    x_ticks = [int(y) for y in x1]
    ax1.set_xticks(x_ticks)
    ax1.set_xticklabels(x_ticks, rotation=45)
    # ax2 for rdd_data
    label2 = 'share of RD&D spending on renewable energy technology'
    ymax2 = ymax1 - y_shift
    ymin2 = ymax2 - t_window
    assert ymin2 >= min(rdd_data['year'].unique()), \
        f"bounds of available data exceeded!"
    rdd_data = rdd_data[rdd_data['year'].isin(range(ymin2, ymax2+1))]
    x2 = rdd_data['year']
    y2 = rdd_data['mean']
    y2err = rdd_data['std']
    ax2.errorbar(
        x2, y2, y2err, fmt='-o', elinewidth=elinewidt,
        capsize=capsize, color='C2', label=label2)
    # ax2 formatting
    x_ticks = [int(y) for y in x2]
    ax2.set_xticks(x_ticks)
    ax2.set_xticklabels(x_ticks, rotation=45)
    # fill error
    if kwargs.get('fill_error', False):
        ax1.fill_between(x1, y1+y1err,y1-y1err, color='C0', alpha=0.1)
        ax2.fill_between(x2, y2+y2err,y2-y2err, color='C2', alpha=0.1)
    # global formatting
    default_title = (
        'Share of RD&D budged invested in renewable energy sources\nvs.'+
        '\nShare of renewables in primary energy\ntime shifted'+
        f'by {y_shift} years'
    )
    ax1.set_xlabel(f'Time for {label1}')
    ax2.set_xlabel(f'Time for {label2}')
    ax1.set_ylabel('Share in %') # global y-axis label
    ax1.set_title(kwargs.get('title', default_title))
    lin1, lab1 = ax1.get_legend_handles_labels()
    lin2, lab2, = ax2.get_legend_handles_labels()
    labels = ['\n'.join(wrap(l, 10)) for l in lab1+lab2]
    ax1.legend(
        lin1 + lin2, labels, bbox_to_anchor=(1.0,1.0),
        loc='upper left', labelspacing=1)

    ax1.grid(axis=kwargs.get('gridaxsis', 'y'))
    ax1.set_axisbelow(True)

