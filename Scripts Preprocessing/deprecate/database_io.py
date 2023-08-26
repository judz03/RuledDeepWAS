'''
This library contains code to load the preprocessed data directly
'''
import pandas as pd


def load_db(database_path:str, nrows:int=None):
    '''
    This little function works only for those databases already preprocessed and stored
    in a .csv file created with pandas.
    '''
    path = database_path
    database = pd.read_csv(path,
                          dtype=str,
                          nrows=nrows)
    return database