from typing import List, Optional
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def txt_to_csv(
    fpath : str, 
    colnames : List,
    outpath : str, 
    return_df : bool =False
    ) -> Optional[pd.DataFrame]:
    """
    Converts .txt data files provided by the International Energy Agency (IEA)
    into .csv files, stores these .csv files and optionally returns a pandas
    data frame containing the data.

    Parameters
    ----------
    fpath : path to .txt file
    colnames: list of column names to use in .csv file
    outpath: path to save .csv file at
    retrun_df: whether or not function returns the .csv file data as a pandas
        data frame

    Returns
    -------
    only returns pandas data frame in case return_df is set to True
    """
    with open(fpath, 'r') as file:
        try:
            data = file.readlines()
            # init list for storing data in proper csv format
            fixed_data = []
            for line in data:
                fixedline = []
                for word in line.split(' '):
                    # only append list items that are not the empty string
                    if word != '':
                        fixedline.append(word)
                    else:
                        continue
                # append proper csv line
                fixed_data.append(','.join(fixedline))
        except IOError:
            file.close()
    # create file if not existent
    if not os.path.exists(outpath):
        Path(outpath).touch()
        with open(outpath, 'a') as file:
            file.writelines(fixed_data)
        # sanity check
        assert len(colnames) == len(fixedline) - 1, \
         'length of colnames must match number of columns in data set'
        # read csv file into temporary data frame for convenience of further
        # formatting
        temp_df =  pd.read_csv(
            outpath,
            names=colnames, # add column names
            sep=',',
            index_col=False
        )
        # replace unituitive missing value reperesentation in original data
        # by np.NaN values
        temp_df.replace(
            to_replace='..',
            value= np.NaN,
            inplace=True
        )
        # save as csv
        temp_df.to_csv(outpath, index=False)
    else:
        temp_df = pd.read_csv(
            outpath,
            sep=',',
            index_col=False
        )
    if return_df:
        return temp_df

def get_missing(data : pd.DataFrame, **kwargs) -> pd.DataFrame:
    """
    Retrns pandas data frame containing the number of missing entries (i.e.
    enteries with NaN values) and the years corresponding to these
    missing entries, per country.
    
    Obviously, less missing entreis are better. Also less missing
    entries for adjacent years are better in terms of ability to interpolate
    missing values.

    Parameters
    ----------
    data : pandas data frame containing IEA RD&D data

    Keyword arguments:
    ymin = lower bound of data range to consider
    ymax = upper bound of data range to consider

    Returns
    -------
    pandas data frame containing info about missing data per country
    """
    years = data['year'].unique()
    ymin = kwargs.get('ymin', min(years))
    ymax = kwargs.get('ymax', max(years))
    c = [] # countries
    n_entry_missing = [] # number of entries missing
    year_missing = [] # years corresponding to missing entries
    for country in data['country'].unique():
        sub_data = data[
            (data['country'] == country) &
            (data['year'] >= ymin) &
            (data['year'] <= ymax)
        ]
        n_missing = sub_data['budget'].isna().sum()
        y_missing = sub_data[sub_data['budget'].isna()]['year'].unique()
        c.append(country)
        n_entry_missing.append(n_missing)
        year_missing.append(y_missing)
    missing_df = pd.DataFrame()
    missing_df['country'] = c
    missing_df['n_entry_missing'] = n_entry_missing
    missing_df['year_missing'] = year_missing
    # sort by number of entries missing
    missing_df.sort_values(by='n_entry_missing', inplace=True)
    missing_df.reset_index(drop=True, inplace=True)
    return missing_df

def interpolate(data : pd.DataFrame) -> None:
    """
    Replace missing RD&D budget vaules of single years by mean of
    budget values of adjacent years

    Parameters
    ----------
    data: pandas data frame containing RD&D budget data 
    """
    for i in data[data['budget'].isna()].index:
        country = data.iloc[i]['country']
        year = data.iloc[i]['year']
        flow = data.iloc[i]['flow']
        adj_data = data[
            (data['country'] == country) &
            (data['flow'] == flow) &
            (data['year'].isin([year-1,year+1]))
        ]
        mean = np.mean(adj_data['budget'])
        data.iloc[i,4] = mean