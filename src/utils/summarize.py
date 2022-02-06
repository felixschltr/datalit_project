from typing import List
import numpy as np
import pandas as pd

def select_years(data, country, **kwargs):
    # get years available in sub data set
    y_available = data['year'].unique()
    # default to min. and max. available year in data set
    ymin = kwargs.get('ymin', min(y_available))
    ymax = kwargs.get('ymax', max(y_available))
    for x, s, f in [(ymin, 'ymin', min), (ymax, 'ymax', max)]:
        if not x in y_available:
            x = f(y_available)
            print(f'specified {s[1:]}. year to consider ({s}) not avialable for'\
                 + f' country {country}, falling back to {s[1:]}. available year '\
                 + f'which is {f(y_available)}.')
    return ymin, ymax
 

def get_diffs(nrg_data : pd.DataFrame, **kwargs) -> pd.DataFrame:
    """
    Caluclates the difference in share from renewable sources between the years
    ymax and ymin for each country. If ymin and ymax are not specified, falls
    back to min. available year and max. available year for each country,
    respectively. Behaves similarly, if specified ymin and ymax are out of
    range for the available data.

    Parameters
    ----------
    nrg_data: data frame containing energy data including data about the share
        from renewable energy sources, per country
    
    Keyword arguments:
    ymin: min. year to consider for calulation of difference; defaults to min.
        available year per country
    ymax: max. year to consider for calulation of difference; defaults to max.
        available year per country
    coi: name of column of interest to caluclate difference on; defaults to
        'share_renewables'
    
    Retruns
    -------
    pandas data frame containing diffrenece in share from renewable sources 
    between ymax and ymin per country
    """
    c = [] # storing countries
    diffs = [] # storing differences
    for country in nrg_data['country'].unique():
        sub_data = nrg_data[nrg_data['country'] == country]
        ymin, ymax = select_years(nrg_data, country, **kwargs)
        # select data between ymin and ymax
        sub_data = sub_data[sub_data['year'].isin(range(ymin, ymax + 1))]
        # sort by year
        sub_data = sub_data.sort_values(by='year',ascending=False)
        # specify column of interest
        coi = kwargs.get('coi', 'share_renewables')
        # compute difference between ymax and ymin
        diff = (
            sub_data.iloc[0][coi] - 
            sub_data.iloc[-1][coi]
        )
        c.append(country)
        diffs.append(diff)
    # create new return data frame
    diff_df = pd.DataFrame()
    diff_df['country'] = c
    diff_df['share_difference'] = diffs
    return diff_df

def get_mean_std(data : pd.DataFrame, selected : List[str], coi : str, **kwargs) -> pd.DataFrame:
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
    # specify min. and max. years
    ymin = kwargs.get('ymin', min(data['year'].unique()))
    ymax = kwargs.get('ymax', max(data['year'].unique()))
    # filter data by selected countries and years
    data = data[
        (data['country'].isin(selected)) &
        (data['year'].isin(range(ymin, ymax+1)))
        ]
    # perform groupby
    groupby_cols = kwargs.get('groupby_cols', ['year'])
    mean_df = (
        data.groupby(groupby_cols).mean()
        .reset_index()
        [groupby_cols + [coi]]
    )
    std_df = (
        data.groupby(groupby_cols).std(ddof=1)
        .reset_index()
        [groupby_cols + [coi]]
    )
    # merge into single data frame
    merged_df = mean_df.merge(std_df, how='inner', on=groupby_cols)
    # rename columns
    merged_df.rename(columns={f'{coi}_x':'mean', f'{coi}_y':'std'}, inplace=True)
    return merged_df