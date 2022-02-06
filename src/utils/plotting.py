import os
import numpy as np
import pandas as pd
import cycler
import matplotlib as mpl
import matplotlib.pyplot as plt
from textwrap import wrap
import tueplots
from tueplots import bundles
from .summarize import get_diffs, get_mean_std

# tueplots settings
plt.rcParams.update(bundles.neurips2021(usetex=True, family="serif"))
default_fs = tueplots.bundles.neurips2021()['figure.figsize']
textsize = tueplots.bundles.neurips2021()['xtick.labelsize']
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
    """
    Produces scatter plot with mean share of R&D speniding on renewable energy
    technology on x-axis and difference in share of primary energy from renewable
    sources on the y-axis. 

    Parameters
    ----------
    rdd_data : pandas data frame containing data on R&D spending
    ymin_rdd : min. year to consider for computation of mean
    ymax_rdd : max. year to consider for computation of mean
    coi_rdd : name of column of interest to compute mean over
    nrg_data : pandas data frame containing data on primary energy consumption
    ymin_nrg : min. year to consider for computation of difference
    ymax_nrg : max. year to consider for computation of difference
    coi_nrg : name of column of interest to compute difference over

    Keyword arguments
    figsize: {tuple of scalar} size of plot
    selected: {list of str} selected countries for plotting; defaults to all
    x_corr_shift : {scalar} shift of correlation annotation along x-axis
    y_corr_shift : {scalar} shift of correlation annotation along y-axis
    tilte : {str} plot title
    y_label : {str} y-axis label
    x_label : {str} x-axis label
    gridaxis : ['y','x','both'] visibility of grid; defaults to 'y'
    save : {bool} whether plot should be saved as .pdf; defaults to False
    path : {str} base directroy for saving the plot

    """
    fig, ax = plt.subplots(figsize=kwargs.get('figsize', default_fs))
    # get differences in share of renewables in energy data
    diff_df = get_diffs(nrg_data, ymin=ymin_nrg, ymax=ymax_nrg, coi=coi_nrg)
    # get means in share of renewables in RD&D dat
    mean_df = get_mean_std(
        rdd_data, rdd_data['country'].unique(),coi=coi_rdd,
        ymin=ymin_rdd, y_max=ymax_rdd, groupby_cols=['country'])
    # merge into singel data frame
    total_df = diff_df.merge(mean_df, how='inner', on='country')
    # optional selection of countries 
    selected = kwargs.get('selected', rdd_data['country'].unique())
    total_df = total_df[total_df['country'].isin(selected)]
    # calculate overall correlation
    r = np.corrcoef(total_df['mean'], total_df['share_difference'])[0][1]
    # acutal plotting
    for country in selected:
        sub_data = total_df[total_df['country'] == country]
        x = sub_data['mean']
        xerr = sub_data['std']
        y = sub_data['share_difference']
        ax.errorbar(x, y, label=country, fmt='-o', elinewidth=1, capsize=2)
        ## add country annotations
        y_offset = (0.95 if country in ['France', 'Norway'] else - 1.3)
        ax.text(
            x - 1.5,
            y + y_offset,
            country,
            fontsize = textsize
        )
     # add correlation annotation
    if kwargs.get('correl', False):
        x_corr = ax.get_xticks()[0]
        y_corr = ax.get_yticks()[-2]
        x_corr_shift = kwargs.get('x_corr_shift', 0)
        y_corr_shift = kwargs.get('y_corr_shift', 2)
        ax.text(
            x_corr + x_corr_shift,
            y_corr + y_corr_shift,
            f'Pearson $r$ = {round(r, 2)}',
            fontsize=9
        )
    # global formatting
    default_title = 'RD\&D spending vs. share of energy from renewable sources'
    ax.set_title(
        '\n'.join(wrap(kwargs.get('title', default_title), 90))
    )
    default_ylab = ('Difference in share of primary energy consumption from ' +
    f'renewable sources between {ymin_nrg} and {ymax_nrg} (in \%)')
    ax.set_ylabel(
        '\n'.join(wrap(kwargs.get('x_label', default_ylab), 60))
    )
    default_xlab = ('Mean share of total energy technology RD\&D ' +
    f'budget invested in renewable sources between {ymin_rdd} and {ymax_rdd}'+ 
    ' (in \%)')
    ax.set_xlabel(
        '\n'.join(wrap(kwargs.get('y_label', default_xlab), 90))
    )
    ax.grid(axis=kwargs.get('gridaxis', 'both'))
    ax.set_axisbelow(True)
    # save settings
    if kwargs.get('save', False):
            # set directory for saving
            dirname = os.path.dirname(__file__)
            save_loc = os.path.join(dirname, kwargs.get('path', '../../gfx'))
            # create if not already existing
            if not os.path.exists(save_loc):
                os.makedirs(save_loc)
            fname = os.path.join(save_loc, 'rdd_vs_nrg_scatter_total.pdf')
            plt.savefig(fname=fname, bbox_inches='tight', dpi=300)

def plot_means_vs_means(
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
    """
        Produces scatter plot with mean share of R&D speniding on renewable energy
        technology on x-axis and difference in share of primary energy from renewable
        sources on the y-axis. 
    
        Parameters
        ----------
        rdd_data : pandas data frame containing data on R&D spending
        ymin_rdd : min. year to consider for computation of mean
        ymax_rdd : max. year to consider for computation of mean
        coi_rdd : name of column of interest to compute mean over
        nrg_data : pandas data frame containing data on primary energy consumption
        ymin_nrg : min. year to consider for computation of mean
        ymax_nrg : max. year to consider for computation of mean
        coi_nrg : name of column of interest to compute mean over
    
        Keyword arguments
        figsize: {tuple of scalar} size of plot
        selected: {list of str} selected countries for plotting; defaults to all
        x_corr_shift : {scalar} shift of correlation annotation along x-axis
        y_corr_shift : {scalar} shift of correlation annotation along y-axis
        tilte : {str} plot title
        y_label : {str} y-axis label
        x_label : {str} x-axis label
        gridaxis : ['y','x','both'] visibility of grid; defaults to 'y'
        save : {bool} whether plot should be saved as .pdf; defaults to False
        path : {str} base directroy for saving the plot
    
        """

    fig, ax = plt.subplots(figsize=kwargs.get('figsize', default_fs))
    # get means in share of renewables in energy data
    nrg_means = get_mean_std(
        nrg_data, nrg_data['country'].unique(),coi=coi_nrg,
        ymin=ymin_nrg, y_max=ymax_nrg, groupby_cols=['country'])
    # get means in share of renewables in RD&D data
    rdd_means = get_mean_std(
        rdd_data, rdd_data['country'].unique(),coi=coi_rdd,
        ymin=ymin_rdd, y_max=ymax_rdd, groupby_cols=['country'])
    # merge into singel data frame
    total_df = nrg_means.merge(rdd_means, how='inner', on='country', suffixes=('_nrg', '_rdd'))
    # optional selection of countries 
    selected = kwargs.get('selected', rdd_data['country'].unique())
    total_df = total_df[total_df['country'].isin(selected)]
    # calculate overall correlation
    r = np.corrcoef(total_df['mean_nrg'], total_df['mean_rdd'])[0][1]
    # acutal plotting
    for country in selected:
        sub_data = total_df[total_df['country'] == country]
        x = sub_data['mean_rdd']
        xerr = sub_data['std_rdd']
        y = sub_data['mean_nrg']
        ax.errorbar(x, y, label=country, fmt='-o', elinewidth=1, capsize=2)
        ## add country annotations
        y_offset = (0.95 if country in ['France', 'Norway'] else - 1.3)
        ax.text(
            x - 1.5,
            y + y_offset,
            country,
            fontsize = textsize
        )
     # add correlation annotation
    if kwargs.get('correl', False):
        x_corr = ax.get_xticks()[0]
        y_corr = ax.get_yticks()[-2]
        x_corr_shift = kwargs.get('x_corr_shift', 0)
        y_corr_shift = kwargs.get('y_corr_shift', 2)
        ax.text(
            x_corr + x_corr_shift,
            y_corr + y_corr_shift,
            f'Pearson $r$ = {round(r, 2)}',
            fontsize=9
        )
     # global formatting
    default_title = 'RD\&D spending vs. share of energy from renewable sources'
    ax.set_title(
        '\n'.join(wrap(kwargs.get('title', default_title), 90))
    )
    default_ylab = ('Mean share of primary energy consumption from renewable' +
    f' sources between {ymin_nrg} and {ymax_nrg} (in \%)')
    ax.set_ylabel(
        '\n'.join(wrap(kwargs.get('x_label', default_ylab), 60))
    )
    default_xlab = ('Mean share of total energy technology RD\&D ' +
    f'budget invested in renewable sources between {ymin_rdd} and {ymax_rdd}'+ 
    ' (in \%)')
    ax.set_xlabel(
        '\n'.join(wrap(kwargs.get('y_label', default_xlab), 90))
    ) 
    ax.grid(axis=kwargs.get('gridaxis', 'both'))
    ax.set_axisbelow(True)
    # saving settings
    if kwargs.get('save', False):
            # set directory for saving
            dirname = os.path.dirname(__file__)
            save_loc = os.path.join(dirname, kwargs.get('path', '../../gfx'))
            # create if not already existing
            if not os.path.exists(save_loc):
                os.makedirs(save_loc)
            fname = os.path.join(save_loc, 'rdd_vs_nrg_means_scatter_total.pdf')
            plt.savefig(fname=fname, bbox_inches='tight', dpi=300)

def plot_diffs_vs_diffs(
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
    """
        Produces scatter plot with difference share of R&D speniding on renewable energy
        technology on x-axis and difference in share of primary energy from renewable
        sources on the y-axis. 
    
        Parameters
        ----------
        rdd_data : pandas data frame containing data on R&D spending
        ymin_rdd : min. year to consider for computation of differnece
        ymax_rdd : max. year to consider for computation of differnece
        coi_rdd : name of column of interest to compute differnece over
        nrg_data : pandas data frame containing data on primary energy consumption
        ymin_nrg : min. year to consider for computation of differnece
        ymax_nrg : max. year to consider for computation of differnece
        coi_nrg : name of column of interest to compute differnece over
    
        Keyword arguments
        figsize: {tuple of scalar} size of plot
        selected: {list of str} selected countries for plotting; defaults to all
        correl : {bool} wheter to plot correlation annotation
        x_corr_shift : {scalar} shift of correlation annotation along x-axis
        y_corr_shift : {scalar} shift of correlation annotation along y-axis
        tilte : {str} plot title
        y_label : {str} y-axis label
        x_label : {str} x-axis label
        gridaxis : ['y','x','both'] visibility of grid; defaults to 'y'
        save : {bool} whether plot should be saved as .pdf; defaults to False
        path : {str} base directroy for saving the plot
    
        """


    fig, ax = plt.subplots(figsize=kwargs.get('figsize', default_fs))
    # get diffs in share of renewables in energy data
    nrg_diffs = get_diffs(nrg_data, ymin=ymin_nrg, ymax=ymax_nrg, coi=coi_nrg)
    # get diffs in share of renewables in RD&D data
    rdd_diffs = get_diffs(rdd_data, ymin=ymin_rdd, ymax=ymax_rdd, coi=coi_rdd)
    # merge into singel data frame
    total_df = nrg_diffs.merge(rdd_diffs, how='inner', on='country', suffixes=('_nrg', '_rdd'))
    # optional selection of countries 
    selected = kwargs.get('selected', rdd_data['country'].unique())
    total_df = total_df[total_df['country'].isin(selected)]
    # calculate overall correlation
    r = np.corrcoef(total_df['share_difference_nrg'], total_df['share_difference_rdd'])[0][1]
    # acutal plotting
    for country in selected:
        sub_data = total_df[total_df['country'] == country]
        x = sub_data['share_difference_rdd']
        y = sub_data['share_difference_nrg']
        ax.errorbar(x, y, label=country, fmt='-o', elinewidth=1, capsize=2)
        ## add country annotations
        y_offset = (0.95 if country in ['France', 'Norway'] else - 1.3)
        ax.text(
            x - 1.5,
            y + y_offset,
            country,
            fontsize = textsize
        )
     # add correlation annotation
    if kwargs.get('correl', False):
        x_corr = ax.get_xticks()[0]
        y_corr = ax.get_yticks()[-2]
        x_corr_shift = kwargs.get('x_corr_shift', 0)
        y_corr_shift = kwargs.get('y_corr_shift', 2)
        ax.text(
            x_corr + x_corr_shift,
            y_corr + y_corr_shift,
            f'Pearson $r$ = {round(r, 2)}',
            fontsize=9
        )
    # global formatting
    default_title = 'RD\&D spending vs. share of energy from renewable sources'
    ax.set_title(
        '\n'.join(wrap(kwargs.get('title', default_title), 90))
    )
    default_ylab = ('Difference in share of primary energy consumption from ' +
    f'renewable sources between {ymin_nrg} and {ymax_nrg} (in \%)')
    ax.set_ylabel(
        '\n'.join(wrap(kwargs.get('x_label', default_ylab), 50))
    )
    default_xlab = ('Difference in share of total energy technology RD\&D ' +
    f'budget invested in renewable sources between {ymin_rdd} and {ymax_rdd}'+ 
    ' (in \%)')
    ax.set_xlabel(
        '\n'.join(wrap(kwargs.get('y_label', default_xlab), 90))
    )
    ax.grid(axis=kwargs.get('gridaxis', 'both'))
    ax.set_axisbelow(True)
    # saving settings
    if kwargs.get('save', False):
            # set directory for saving
            dirname = os.path.dirname(__file__)
            save_loc = os.path.join(dirname, kwargs.get('path', '../../gfx'))
            # create if not already existing
            if not os.path.exists(save_loc):
                os.makedirs(save_loc)
            fname = os.path.join(save_loc, 'rdd_vs_nrg_means_scatter_total.pdf')
            plt.savefig(fname=fname, bbox_inches='tight', dpi=300)



def plot_rdd(rdd_data : pd.DataFrame, **kwargs):
    """
    Create plot of RD\&D data for each country with year on the x-axis and
    share of budget invested in renewable energy sources on y-axis.

    Parameters
    ----------
    rdd_data: data on RD\&D spending per country

    Keyword arguments:
    avg_window: {int} year window for computing rolling average
        defaults to 1, i.e. no averaging
    figsize: {tuple of int} defaults to tuebplots neurips21 figsize
    selected: {list/iterable of strings} selection of countries to plot
    plot_non_avg {boolean}: when True, plots true signal in background
    alpha {scalar} : transparency of ture signal (if plotted)
    gridaxis {str} : defaults to 'y'

    """
    fig, ax = plt.subplots(figsize=kwargs.get('figsize', default_fs))
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
    default_title = ('Share of energy technology RD\&D budget invested in' +
    'renewable energy sources by country and  year')
    ax.set_title(
        '\n'.join(wrap(kwargs.get('title', default_title), 60))
    )
    ax.set_xlabel('Year')
    ylab = 'Share of RD\&D budget invested in renewable energy sources (in \%)'
    if avg_window != 1:
        ylab += f' [{avg_window}-year running average]'
    ax.set_ylabel('\n'.join(wrap(ylab, 50)))
    ax.legend(bbox_to_anchor=(1.0,1.0), loc='upper left')
    ax.grid(axis=kwargs.get('gridaxis', 'y'))
    ax.set_axisbelow(True)

def plot_lines(
    rdd_data : pd.DataFrame,
    nrg_data : pd.DataFrame,
    t_window : int = 10,
    y_shift : int = 10,
    avg_window : int = 1,
    **kwargs
    ):
    """
    Plots data of R&D spending on renewable energy technology and data on
    primary energy consumption from renewable sources as separete lines. If
    time shift is specified, sequence of R&D data is shifted forward relative
    to sequence of primary energy data and two separate x-axes are plotted.

    Parameters
    ----------
    rdd_data : pandas data frame containing data on R&D spending
    nrg_data : pandas data frame containing data on primary energy consumption
    t_window : time window of data to consider (in year) going backward from
        latest availabe year of energy data
    y_shift : foward shift of R&D data relative to energy datd
    avg_window: year window for computing rolling average
        defaults to 1, i.e. no averaging

    Keyword arguments:
    figsize: {tuple of scalar} size of plot
    selected: {list of str} selected countries for plotting; defaults to all
    ymax1 : {int} maximum date to consider for energy data
    plot_non_avg : {bool} whether to plot non-averaged data; defaults to False
    alpha : {scalar} transparency of non-averaged signal; defaults to 0.6
    correl : {bool} whether to plot correlation coefficient; defaults to False
    x_corr_shift : {scalar} shift of correlation annotation along x-axis
    y_corr_shift : {scalar} shift of correlation annotation along y-axis
    tilte : {str} plot title
    gridaxis : ['y','x','both'] visibility of grid; defaults to 'y'
    save : {bool} whether plot should be saved as .pdf; defaults to False
    path : {str} base directroy for saving the plot
   """
    fig, ax1 = plt.subplots(figsize=kwargs.get('figsize', default_fs))
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
    selected = kwargs.get('selected', rdd_data['country'].unique())
    for country in selected:
        # plot nrg_data
        label1_base = 'share of primary energy consumption from renewable sources'
        sub_data = (
            nrg_data[(nrg_data['country'] == country)]
            .sort_values(by='year')
        )
        x1 = sub_data['year']
        y1 = sub_data['prim_nrg_share_renewables']
        ax1.plot(x1, y1, '-s', c='C0', label=label1_base)
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
        label2_base = 'share of RD\&D spending on renewable energy technology'
        if avg_window != 1:
            label2 = label2_base + f' ({avg_window}-year rolling average)'
        if kwargs.get('plot_non_avg', False):
            ax2.plot(x2, y2, ls='--',alpha=kwargs.get('alpha', 0.6), c='C2', label='share RD\&D non-averaged')
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
            x_corr = ax1.get_xticks()[0]
            y_corr = ax1.get_yticks()[-2]
            x_corr_shift = kwargs.get('x_corr_shift', 0)
            y_corr_shift = kwargs.get('y_corr_shift', 2)
            ax1.text(
                x_corr + x_corr_shift,
                y_corr - y_corr_shift,
                f'Pearson $r$ = {round(pearson_r, 3)}',
                fontsize = 9
                )
    # global formatting
    default_title = (
        'Share of RD\&D budged invested in renewable energy sources vs. '+
        'share of renewables in primary energy consumption time shifted '+
        f'by {y_shift} years'
    )
    ax1.set_xlabel(f'Time for {label1_base}')
    ax2.set_xlabel(f'Time for {label2_base}')
    ax1.set_ylabel('Share in \%') # global y-axis label
    ax1.set_title(
        '\n'.join(wrap(kwargs.get('title', default_title), 90))
    )
    # create single legend
    lin1, lab1 = ax1.get_legend_handles_labels()
    lin2, lab2, = ax2.get_legend_handles_labels()
    labels = ['\n'.join(wrap(l, 10)) for l in lab1+lab2]
    ax1.legend(
        lin1 + lin2, labels, bbox_to_anchor=(1.0,1.0),
        loc='upper left', labelspacing=1)
    ax1.grid(axis=kwargs.get('gridaxsis', 'y'))
    ax1.set_axisbelow(True)
    if kwargs.get('save', False):
        # set directory for saving
        dirname = os.path.dirname(__file__)
        save_loc = os.path.join(dirname, kwargs.get('path', '../../gfx'))
        # create if not already existing
        if not os.path.exists(save_loc):
            os.makedirs(save_loc)
        name = '_'.join(selected)
        fname = os.path.join(save_loc, f'rdd_vs_nrg_lines_{name}.pdf')
        plt.savefig(fname=fname, bbox_inches='tight', dpi=300)

def plot_lines_total(
    rdd_data : pd.DataFrame,
    nrg_data : pd.DataFrame,
    t_window : int = 10,
    y_shift : int = 10,
    **kwargs
    ):
    """
    Plots data of mean R&D spending on renewable energy technology and data on
    mean primary energy consumption from renewable sources, averaged over 
    countries as separete lines. If time shift is specified, sequence of R&D
    data is shifted forward relative to sequence of primary energy data and
    two separate x-axes are plotted.

    Parameters
    ----------
    rdd_data : pandas data frame containing data on R&D spending
    nrg_data : pandas data frame containing data on primary energy consumption
    t_window : time window of data to consider (in year) going backward from
        latest availabe year of energy data
    y_shift : foward shift of R&D data relative to energy datd
    avg_window: year window for computing rolling average
        defaults to 1, i.e. no averaging

    Keyword arguments:
    elinewidth : {scalar} width of error bar line
    capsize : {scalar} size of error bar caps
    figsize: {tuple of scalar} size of plot
    selected: {list of str} selected countries for plotting; defaults to all
    ymax1 : {int} maximum date to consider for energy data
    plot_non_avg : {bool} whether to plot non-averaged data; defaults to False
    alpha : {scalar} transparency of non-averaged signal; defaults to 0.6
    correl : {bool} whether to plot correlation coefficient; defaults to False
    x_corr_shift : {scalar} shift of correlation annotation along x-axis
    y_corr_shift : {scalar} shift of correlation annotation along y-axis
    tilte : {str} plot title
    gridaxis : ['y','x','both'] visibility of grid; defaults to 'y'
    save : {bool} whether plot should be saved as .pdf; defaults to False
    path : {str} base directroy for saving the plot
   """
 
    assert t_window >= 0, 't_window must be positive!'
    assert y_shift >= 0, 'y_shift must be positive!'
    # formatting of error bars
    elinewidth = kwargs.get('elinewidth', 1)
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
    fig, ax1 = plt.subplots(figsize=kwargs.get('figsize', default_fs))
    # use on x-axis of two depending on parameter y_shift
    ax2 = (ax1 if y_shift == 0 else ax1.twiny()) # create twin y axis
    # ax1 for nrg_data
    label1 = 'share of primary energy consumption from renewable sources'
    ymax1 = kwargs.get('ymax1',max(nrg_data['year'].unique()))
    ymin1 = ymax1 - t_window
    assert ymin1 >= min(nrg_data['year'].unique()), \
        f"bounds of available data exceeded!"
    nrg_data = nrg_data[nrg_data['year'].isin(range(ymin1, ymax1+1))]
    x1 = nrg_data['year']
    y1 = nrg_data['mean']
    y1err = nrg_data['std']
    ax1.errorbar(
        x1, y1, y1err, fmt='-s', elinewidth=elinewidth,
        capsize=capsize, color='C0', label=label1)
    # ax1 formatting
    x_ticks = [int(y) for y in x1]
    ax1.set_xticks(x_ticks)
    ax1.set_xticklabels(x_ticks, rotation=45)
    # ax2 for rdd_data
    label2 = 'share of RD\&D spending on renewable energy technology'
    ymax2 = ymax1 - y_shift
    ymin2 = ymax2 - t_window
    assert ymin2 >= min(rdd_data['year'].unique()), \
        f"bounds of available data exceeded!"
    rdd_data = rdd_data[rdd_data['year'].isin(range(ymin2, ymax2+1))]
    x2 = rdd_data['year']
    y2 = rdd_data['mean']
    y2err = rdd_data['std']
    ax2.errorbar(
        x2, y2, y2err, fmt='-o', elinewidth=elinewidth,
        capsize=capsize, color='C2', label=label2)
    # ax2 formatting
    x_ticks = [int(y) for y in x2]
    ax2.set_xticks(x_ticks)
    ax2.set_xticklabels(x_ticks, rotation=45)
    # fill error
    if kwargs.get('fill_error', False):
        ax1.fill_between(x1, y1+y1err,y1-y1err, color='C0', alpha=0.1)
        ax2.fill_between(x2, y2+y2err,y2-y2err, color='C2', alpha=0.1)
    # add correlation annotation
    if kwargs.get('correl', False):
        r = np.corrcoef(y1, y2)[0][1]
        x_corr = ax1.get_xticks()[0]
        y_corr = ax1.get_yticks()[-2]
        x_corr_shift = kwargs.get('x_corr_shift', 0)
        y_corr_shift = kwargs.get('y_corr_shift', 2)
        ax1.text(
            x_corr + x_corr_shift,
            y_corr + y_corr_shift,
            f'Pearson $r$ = {round(r, 2)}',
            fontsize=9
        )
    # optional formatting
    if y_shift == 0: # applies if signals are not time-shifted
        ax1.set_xlabel('Time')
        lin1, lab1 = ax1.get_legend_handles_labels()
        labels = ['\n'.join(wrap(l, 10)) for l in lab1]
        ax1.legend(
            lin1, labels, bbox_to_anchor=(1.0,1.0),
            loc='upper left', labelspacing=1)
    else: # applies if singals are time-shifted
        ax1.set_xlabel(f'Time for {label1}')
        ax2.set_xlabel(f'Time for {label2}')
        lin1, lab1 = ax1.get_legend_handles_labels()
        lin2, lab2, = ax2.get_legend_handles_labels()
        labels = ['\n'.join(wrap(l, 10)) for l in lab1+lab2]
        ax1.legend(
            lin1 + lin2, labels, bbox_to_anchor=(1.0,1.0),
            loc='upper left', labelspacing=1)
    # global formatting
    default_title = (
        'Share of RD\&D budged invested in renewable energy sources vs. '+
        'share of renewables in primary energy consumption time shifted '+
        f'by {y_shift} years'
    )
    ax1.set_title(
        '\n'.join(wrap(kwargs.get('title', default_title), 90))
    )
    ax1.set_ylabel('Share in \%') # global y-axis label
    ax1.grid(axis=kwargs.get('gridaxsis', 'y'))
    ax1.set_axisbelow(True)
    if kwargs.get('save', False):
        # set directory for saving
        dirname = os.path.dirname(__file__)
        save_loc = os.path.join(dirname, kwargs.get('path', '../../gfx'))
        # create if not already existing
        if not os.path.exists(save_loc):
            os.makedirs(save_loc)
        name = f'rdd_vs_nrg_lines_total_twindow_{t_window}_yshift_{y_shift}.pdf'
        fname = os.path.join(save_loc, name)
        plt.savefig(fname=fname, bbox_inches='tight', dpi=300)