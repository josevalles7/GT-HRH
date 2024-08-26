"""
Functions for HydroSOS Status product
Based on STATUSCALCV2.R by Katie Facer-Childs and Ezra Kitson
@author: Jose Valles (09/11/2023)
DINAGUA - URUGUAY
"""
import pandas as pd
import numpy as np
import os

# ------------------------------------------------------------------------------
# Hydrological Status Parameters
# ------------------------------------------------------------------------------

stdStart = 1981
stdEnd = 2010
percentile = [0.10, 0.25, 0.75, 0.90]
max_pct_missing = 50

values = ['High flow','Above normal','Normal range','Below normal','Low flow']
flow_cat = [5,4,3,2,1]

# ------------------------------------------------------------------------------
# Main functions
# ------------------------------------------------------------------------------
   
def monthly_status(DISCHARGE_MONTHLY):
    # Calculate long-term average
    DISCHARGE_SELECTION = DISCHARGE_MONTHLY[(DISCHARGE_MONTHLY['year'] >= stdStart) & (DISCHARGE_MONTHLY['year'] <= stdEnd)]
    DISCHARGE_AVERAGE = DISCHARGE_SELECTION.groupby(DISCHARGE_SELECTION.month).mean()
    DISCHARGE_AVERAGE = DISCHARGE_AVERAGE.reindex(columns=['mean_flow'])
    # Calcular indicadores
    DISCHARGE_MONTHLY['percentile_flow'] = np.nan
    DISCHARGE_MONTHLY['rank_average'] = np.nan
    DISCHARGE_MONTHLY['complete%'] = np.nan

    for i in range(len(DISCHARGE_MONTHLY)):
        # Extract the current month 
        m = DISCHARGE_MONTHLY.month[i]
        # Extract the current year
        y = DISCHARGE_MONTHLY.year[i]
        DISCHARGE_MONTHLY.loc[DISCHARGE_MONTHLY.eval('month==@m & year==@y'),'rank_average']  = DISCHARGE_MONTHLY.query('month==@m')['mean_flow'].rank()
        DISCHARGE_MONTHLY.loc[DISCHARGE_MONTHLY.eval('month==@m & year==@y'),'complete%']  = DISCHARGE_MONTHLY.query('month==@m')["mean_flow"].notnull().sum()
        DISCHARGE_MONTHLY.loc[DISCHARGE_MONTHLY.eval('month==@m & year==@y'),'percentile_flow'] = (DISCHARGE_MONTHLY['mean_flow'][i] - DISCHARGE_AVERAGE.query('month == @m')["mean_flow"].item()) / DISCHARGE_AVERAGE.query('month == @m')["mean_flow"].item()

    DISCHARGE_MONTHLY['weibell_rank'] = DISCHARGE_MONTHLY['rank_average']/(DISCHARGE_MONTHLY['complete%']+1)

    criteria = [DISCHARGE_MONTHLY['weibell_rank'].between(percentile[3],1.00),
        DISCHARGE_MONTHLY['weibell_rank'].between(percentile[2],percentile[3]),
        DISCHARGE_MONTHLY['weibell_rank'].between(percentile[1],percentile[2]),
        DISCHARGE_MONTHLY['weibell_rank'].between(percentile[0],percentile[1]),
        DISCHARGE_MONTHLY['weibell_rank'].between(0.00,percentile[0])]

    DISCHARGE_MONTHLY['percentile_range'] = np.select(criteria,values,None)
    DISCHARGE_MONTHLY['flowcat'] = np.select(criteria,flow_cat,pd.NA)

    return DISCHARGE_MONTHLY

def quarterly_status(DISCHARGE_THREE_MONTHS):
    # # Calculate long-term average
    DISCHARGE_SELECTION_THREE_MONTH = DISCHARGE_THREE_MONTHS[(DISCHARGE_THREE_MONTHS['year'] >= stdStart) & (DISCHARGE_THREE_MONTHS['year'] < stdEnd)]
    DISCHARGE_AVERAGE_THREE_MONTH = DISCHARGE_SELECTION_THREE_MONTH.groupby(DISCHARGE_SELECTION_THREE_MONTH.startMonth).mean()
    DISCHARGE_AVERAGE_THREE_MONTH = DISCHARGE_AVERAGE_THREE_MONTH.reindex(columns=['mean_flow'])
    # Calcular indicadores
    DISCHARGE_THREE_MONTHS['percentage_flow'] = np.nan
    DISCHARGE_THREE_MONTHS['rank_average'] = np.nan
    DISCHARGE_THREE_MONTHS['complete%'] = np.nan

    for i in range(len(DISCHARGE_THREE_MONTHS)):
        # Extract the current month 
        m = DISCHARGE_THREE_MONTHS.startMonth[i]
        # Extract the current year
        y = DISCHARGE_THREE_MONTHS.year[i]
        DISCHARGE_THREE_MONTHS.loc[DISCHARGE_THREE_MONTHS.eval('startMonth==@m & year==@y'),'rank_average']  = DISCHARGE_THREE_MONTHS.query('startMonth==@m')['mean_flow'].rank()
        DISCHARGE_THREE_MONTHS.loc[DISCHARGE_THREE_MONTHS.eval('startMonth==@m & year==@y'),'complete%']  = DISCHARGE_THREE_MONTHS.query('startMonth==@m')["mean_flow"].notnull().sum()
        DISCHARGE_THREE_MONTHS.loc[DISCHARGE_THREE_MONTHS.eval('startMonth==@m & year==@y'),'percentage_flow'] = (DISCHARGE_THREE_MONTHS['mean_flow'][i] - DISCHARGE_AVERAGE_THREE_MONTH.query('startMonth == @m')["mean_flow"].item()) / DISCHARGE_AVERAGE_THREE_MONTH.query('startMonth == @m')["mean_flow"].item()

    DISCHARGE_THREE_MONTHS['weibell_rank'] = DISCHARGE_THREE_MONTHS['rank_average']/(DISCHARGE_THREE_MONTHS['complete%']+1)
   
    criteria_three_months = [DISCHARGE_THREE_MONTHS['weibell_rank'].between(percentile[3],1.00),
        DISCHARGE_THREE_MONTHS['weibell_rank'].between(percentile[2],percentile[3]),
        DISCHARGE_THREE_MONTHS['weibell_rank'].between(percentile[1],percentile[2]),
        DISCHARGE_THREE_MONTHS['weibell_rank'].between(percentile[0],percentile[1]),
        DISCHARGE_THREE_MONTHS['weibell_rank'].between(0.00,percentile[0])]

    DISCHARGE_THREE_MONTHS['percentile_range'] = np.select(criteria_three_months,values,None)
    DISCHARGE_THREE_MONTHS['flowcat'] = np.select(criteria_three_months,flow_cat,pd.NA)

    row_labels = {1:'JFM',
            2:'FMA',
            3:'MAM',
            4:'AMJ',
            5:'MJJ',
            6:'JJA',
            7:'JAS',
            8:'ASO',
            9:'SON',
            10:'OND',
            11:'NDE',
            12:'DEF'}
    
    DISCHARGE_THREE_MONTHS['period'] = DISCHARGE_THREE_MONTHS['startMonth'].replace(row_labels) 
    return DISCHARGE_THREE_MONTHS

def annualy_status(DISCHARGE_TWELVE_MONTHS):
    DISCHARGE_SELECTION_TWELVE_MONTH = DISCHARGE_TWELVE_MONTHS[(DISCHARGE_TWELVE_MONTHS['year'] >= stdStart) & (DISCHARGE_TWELVE_MONTHS['year'] < stdEnd)]
    DISCHARGE_AVERAGE_TWELVE_MONTH = DISCHARGE_SELECTION_TWELVE_MONTH.groupby(DISCHARGE_SELECTION_TWELVE_MONTH.startMonth).mean()
    DISCHARGE_AVERAGE_TWELVE_MONTH = DISCHARGE_AVERAGE_TWELVE_MONTH.reindex(columns=['mean_flow'])
    # Calcular indice
    DISCHARGE_TWELVE_MONTHS['percentage_flow'] = np.nan
    DISCHARGE_TWELVE_MONTHS['rank_average'] = np.nan
    DISCHARGE_TWELVE_MONTHS['complete%'] = np.nan

    for i in range(len(DISCHARGE_TWELVE_MONTHS)):
        # Extract the current month 
        m = DISCHARGE_TWELVE_MONTHS.startMonth[i]
        # Extract the current year
        y = DISCHARGE_TWELVE_MONTHS.year[i]
        DISCHARGE_TWELVE_MONTHS.loc[DISCHARGE_TWELVE_MONTHS.eval('startMonth==@m & year==@y'),'rank_average']  = DISCHARGE_TWELVE_MONTHS.query('startMonth==@m')['mean_flow'].rank()
        DISCHARGE_TWELVE_MONTHS.loc[DISCHARGE_TWELVE_MONTHS.eval('startMonth==@m & year==@y'),'complete%']  = DISCHARGE_TWELVE_MONTHS.query('startMonth==@m')["mean_flow"].notnull().sum()
        DISCHARGE_TWELVE_MONTHS.loc[DISCHARGE_TWELVE_MONTHS.eval('startMonth==@m & year==@y'),'percentage_flow'] = (DISCHARGE_TWELVE_MONTHS['mean_flow'][i] - DISCHARGE_AVERAGE_TWELVE_MONTH.query('startMonth == @m')["mean_flow"].item()) / DISCHARGE_AVERAGE_TWELVE_MONTH.query('startMonth == @m')["mean_flow"].item()

    DISCHARGE_TWELVE_MONTHS['weibull_rank'] = DISCHARGE_TWELVE_MONTHS['rank_average']/(DISCHARGE_TWELVE_MONTHS['complete%']+1)
    
    criteria_twelve_months = [DISCHARGE_TWELVE_MONTHS['weibull_rank'].between(percentile[3], 1.00),
                          DISCHARGE_TWELVE_MONTHS['weibull_rank'].between(percentile[2], percentile[3]),
                          DISCHARGE_TWELVE_MONTHS['weibull_rank'].between(percentile[1], percentile[2]),
                          DISCHARGE_TWELVE_MONTHS['weibull_rank'].between(percentile[0], percentile[1]),
                          DISCHARGE_TWELVE_MONTHS['weibull_rank'].between(0.00, percentile[0])]
        
    DISCHARGE_TWELVE_MONTHS['percentile_range'] = np.select(criteria_twelve_months,values,None)
    DISCHARGE_TWELVE_MONTHS['flowcat'] = np.select(criteria_twelve_months,flow_cat,pd.NA)

    return DISCHARGE_TWELVE_MONTHS

def export_csv(groupBy,output_directory,filename):
    groupBy['date'] = pd.to_datetime(groupBy[['year', 'month']].assign(DAY=1))
    groupBy['date'] = groupBy['date'].dt.strftime('%Y-%m-%d')
    groupBy['flowcat'] = groupBy['flowcat'].astype('Int64')
    groupBy.sort_values(['year','month']).filter(['date','flowcat']).to_csv(f"{output_directory}cat_{filename}.csv", index=False)
    print('CSV FILE GENERATED')

def csv_to_json(input_csv,outpu_json):
    allFilesDF = pd.DataFrame()
    # read the CSV files in the data directory
    for index, filename in enumerate(os.listdir(input_csv)):
            with open(input_csv+'/'+filename, mode="r") as fr:
                if filename.endswith('.csv'):
                    df = pd.read_csv(fr)
                    filename = os.path.splitext(str(filename))[0] #remove file extenstion
                    stationID = filename.split('_')[1] #remove cat_
                    df['stationID'] = stationID
                    allFilesDF = pd.concat([allFilesDF,df])

    allFilesDF['date'] = pd.to_datetime(allFilesDF['date'])
    allFilesDF = allFilesDF.sort_values(by='date', ascending=True)
    allFilesDF.drop_duplicates(inplace=True)
    allFilesDF.set_index(['date'], inplace=True)
    allFilesDF.rename(columns={"flowcat":"category"}, inplace=True)
    allFilesDF['category'] = allFilesDF['category'].astype('Int64')

    for date in allFilesDF.index:
        #this happens if there is only one record
        if type(allFilesDF.loc[date]) == pd.core.series.Series:
            pd.DataFrame(allFilesDF.loc[date]).T.to_json(f"{outpu_json}/{date.strftime('%Y-%m')}.json", orient = 'records')
        else: 
            allFilesDF.loc[date].to_json(f"{outpu_json}/{date.strftime('%Y-%m')}.json", orient = 'records')

    print('JSON FILE GENERATED')