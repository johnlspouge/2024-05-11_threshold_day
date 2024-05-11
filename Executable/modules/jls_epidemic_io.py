#!/usr/bin/env python
"""
Contains auxiliary subroutines for user epidemic input, with extensive checking. 
"""
from os.path import isfile
import numpy as np
import pandas as pd
import csv
from io import StringIO
from math import isclose

from jls_callback import Callback
from jls_epidemic_stochastic import Epidemic
from jls_epidemic_compartment import Compartment
from jls_epidemic_watcher_threshold import to_passive_watchers, active_for_difference, active_for_total

class Io:
    DECAY_PARAMS = ['mean','dispersion','prob']
    INFECT_PARAMS = ['R_0']
    CSV = '.csv'
    # Returns a list of dfs with same columns as df but with one row only from df.
    #   The list is useful in itertools to get all combinations of initial counts, etc.
    def to_dfs(df):
        compartment_names = list(df.columns)
        dfs = [None]*len(df)
        for i in range(len(df)):
            counts = df.iloc[i,:].values.tolist()
            dfs[i] = pd.DataFrame([counts], columns=compartment_names)
        return dfs # Returns dfs with the original order of the rows.
    # Returns a dataframe whose column titles are compartments and whose rows are initial counts.
    #   Converts CSV file with column[0]=compartments and columns=sets of initial counts to dataframe.
    #   Compartment names must contain only letters & underscores.
    def dfs_initialize_from(
            ifn, # input CSV filename for dataframe 
            sep=','): # separator for CSV dataframe
        # Checks the contents of a CSV file.
        if not isfile( f'{ifn}' ):
            raise ValueError(f'Error: a valid INPUT_FILE_CSV "{ifn}" is required.')
        try:
            with open(ifn, 'r') as ifh:
                string = ifh.read()
        except:
            raise ValueError( f'Error: Read of file {ifn} failed.' )
        # Adjusts the contents of the dataframe.
        dfs = Io.to_dfs_initialize_from(string, sep=sep) # Permits multiple initial conditions.
        return dfs 
    # Returns a list of dataframes whose column titles are compartments and whose rows are initial counts.
    #   Converts CSV file with column[0]=compartments and columns=sets of initial counts to dataframe.
    #   Compartment names must contain only letters & underscores.
    def to_dfs_initialize_from(
            string, # contents, e.g., of input CSV filename for dataframe 
            sep=','): # separator for CSV dataframe
        # Checks the contents of a CSV file.
        try:
            ifh = StringIO(string)
            reader = csv.reader(ifh)
            rows = []
            for i,row in enumerate(reader):
                if i == 0:
                    continue
                rows.append(row)
                if len(row) != len(rows[0]):
                    raise ValueError(f'Error: All rows must have equal length - "{row}".')
            if not rows:
                raise ValueError(f'Error : empty csv string')
            transposes = [ [] for cell in rows[0]]
            for j,row in enumerate(rows):
                for i,cell in enumerate(row):
                    if j > 0 and cell == '': # A count for a compartment (j==0) is blank.
                        cell = '0'
                    transposes[i].append(cell)
            transposes = list(map(lambda x: ','.join(x), transposes))
            transpose_c = '\n'.join(transposes)
            transpose_c += '\n'
            transpose_io = StringIO(transpose_c)
            df = pd.read_csv(transpose_io, sep=',', dtype=str)
            for i in range(len(df)):
                df.iloc[i] = df.iloc[i].astype(int)
        except:
            raise ValueError( f'Error: Some data has an incorrect type.' )
        # Adjusts the contents of the dataframe.
        Io.check_initialize_from(df)
        dfs = []
        for i in range(len(df)):
            df0 = df.iloc[df.index.isin([i])]
            df0.reset_index(drop=True)
            dfs.append(df0)
        return dfs # a list of dfs with each row of df
    # Checks the DataFrame returned by dfs_initialize_from. 
    def check_initialize_from(df): # dataframe with columns=compartments 1st row=initial counts
        # Checks the contents of the dataframe.
        for i in range(len(df)):
            row = df.loc[i,:]
            lst = row.values.flatten().tolist()
            if row.lt(0).any():
                 raise ValueError(f'Dataframe row {lst} contains a value < 0')
            if len(lst) != len(df.columns):
                raise ValueError('Some compartment names are duplicated.')
    # Converts input file to dataframe with the decays for unstable compartments.
    def df_add_decay_to(
            compartments, # Are decays for valid compartments ?
            ifn, # input filename
            sep=','): # separator for CSV dataframe
        # Checks the contents of a CSV file.
        if not isfile( f'{ifn}' ):
            raise ValueError( f'Error: a valid INPUT_FILE_CSV "{ifn}" is required.' )
        try:
            with open(ifn, 'r') as ifh:
                string = ifh.read()
        except:
            raise ValueError( f'Error: Read of file {ifn} failed.' )
        df = Io.to_df_add_decay_to(compartments, string, sep=sep)
        return df
    # Converts input file to dataframe with the decays for unstable compartments.
    def to_df_add_decay_to(
            compartments, # Are decays for valid compartments ?
            string, # contents, e.g., of input CSV filename for dataframe 
            sep=','): # separator for CSV dataframe
        try:
            ifh = StringIO(string)
            df = pd.read_csv(ifh, sep=sep, dtype={'from':str,'to':str,'mean':str,'dispersion':str,'prob':str})
            # Adjusts the contents of the dataframe.
            df['prob'] = df['prob'].fillna(1.0) # if the CSV cell for 'prob' = '', then the decay distribution is a deterministic atom.
            for c in list(set(df['from'])):
                # mean
                set_of_means = set(df.loc[df['from'] == c]['mean'])-{np.nan}
                if len(set_of_means) == 0:
                    raise ValueError(f'At least one mean decay time for compartment {c} must be a number.')
                elif len(set_of_means) > 1:
                    raise ValueError(f'Every mean decay time for compartment {c} must be equal or blank.')
                (mean,) = set_of_means
                df.loc[(df['from'] == c) & (df['mean'].isna()),'mean'] = mean
                # dispersion
                set_of_dispersions = set(df.loc[df['from'] == c]['dispersion'])-{np.nan}
                if len(set_of_dispersions) > 1:
                    raise ValueError(f'Every dispersion of decay time for compartment {c} must be equal or blank.')
                elif len(set_of_dispersions) == 1:
                    (dispersion,) = set_of_dispersions
                    df.loc[(df['from'] == c) & (df['dispersion'].isna()),'dispersion'] = dispersion
            df['mean'] = df['mean'].astype(float)
            df['dispersion'] = df['dispersion'].astype(float)
            df['prob'] = df['prob'].astype(float)
        except:
            raise ValueError( f'Error: Some data has an incorrect type.' )
        Io.check_add_decay_to(compartments, df)
        return df
    # Adds the decays for unstable compartments.
    def check_add_decay_to(compartments, df): # dataframe with columns ['from','to','mean','dispersion','prob']
        # Checks the contents of the dataframe.
        for c in df['from']:
            if c not in compartments:
                raise ValueError( f'Error: "{c}" in column "from" is not a valid compartment.' )
        for c in df['to']:
            if c not in compartments:
                raise ValueError( f'Error: "{c}" in column "to" is not a valid compartment.' )
        for i, row in df.iterrows():
            if not 0.0 < df.at[i,'mean']:
                raise ValueError(f'Compartment decay mean in row {i} is not 0.0 < {row["mean"]}.')
            if not pd.isnull(df.at[i,'dispersion']) and not 0.0 < df.at[i,'dispersion']:
                raise ValueError(f'Compartment decay mean in row {i} not 0.0 < {row["dispersion"]}.')
            if not 0.0 <= df.at[i,'prob'] <= 1.0:
                raise ValueError(f'Compartment decay probability in row {i} not 0.0 <= {row["prob"]} <= 1.0')
        for c in df['from']:
            if not len(set(df.loc[df['from'] == c]['mean'])) == 1:
                raise ValueError(f'Every mean decay time for compartment {c} must be equal.')
            if len(set(df.loc[df['from'] == c]['dispersion'])) > 1:
                raise ValueError(f'Every dispersion of decay time for compartment {c} must be equal.')
            if not isclose(df.loc[df['from'] == c]['prob'].sum(), 1.0):
                raise ValueError(f'Compartment decay probabilities for {c} do not sum to 1.0.')            
    # Adds the infects for infecting compartments.
    def df_add_infect_to(
            compartments, # Are infections for valid compartments ?
            decay_froms, # Do the infecting compartments decay spontaneously ?
            ifn, # input filename
            sep=','): # separator for CSV dataframe
        # Checks the contents of a CSV file.
        if not isfile( f'{ifn}' ):
            raise ValueError( f'Error: a valid INPUT_FILE_CSV "{ifn}" is required.' )
        try:
            with open(ifn, 'r') as ifh:
                string = ifh.read()
        except:
            raise ValueError( f'Error: Read of file {ifn} failed.' )
        df = Io.to_df_add_infect_to(compartments, decay_froms, string, sep=sep)
        return df
     # Adds the infects for infecting compartments.
    def to_df_add_infect_to(
            compartments, # Are infections for valid compartments ?
            decay_froms, # Do the infecting compartments decay spontaneously ?
            string, # contents, e.g., of input CSV filename for dataframe
            sep=','): # separator for CSV dataframe
        # Checks the contents of a CSV file.
        try:
            ifh = StringIO(string)
            df = pd.read_csv(ifh, sep=sep, dtype={'from':str,'to':str,'infecting':str,'R_0':str})
            df['R_0'] = df['R_0'].astype(float)
            if not len(df):
                raise ValueError( f'Error: There are no infecting states.' )
        except:
            raise ValueError( f'Error: Some data has an incorrect type.' )
        Io.check_add_infect_to(compartments, decay_froms, df)
        return df
    # Checks the infects for infecting compartments.
    def check_add_infect_to(compartments, decay_froms, df): # dataframe with columns ['from','to','infecting','R_0']
        # Checks the contents of the dataframe.
        for c in df['from']:
            if c not in compartments:
                raise ValueError( f'Error: "{c}" in column "from" is not a valid compartment.' )
        for c in df['to']:
            if c not in compartments:
                raise ValueError( f'Error: "{c}" in column "to" is not a valid compartment.' )
        for c in df['infecting']:
            if c not in decay_froms:
                raise ValueError( f'Error: "{c}" in column "infecting" does not decay spontaneously.' )
        for r0 in df['R_0']:
            if not r0 > 0.0:
                raise ValueError( f'Error: "{r0}" in column "R_0" is not positive.' )
    # Adds the infects for infecting compartments.
    def df_watchers(
            compartments, # Are infections for valid compartments ? 
            ifn, # input filename
            sep=','): # separator for CSV dataframe
        # Checks the contents of a CSV file.
        if not isfile( f'{ifn}' ):
            raise ValueError( f'Error: a valid INPUT_FILE_CSV "{ifn}" is required.' )
        try:
            with open(ifn, 'r') as ifh:
                string = ifh.read()
        except:
            raise ValueError( f'Error: Some data has an incorrect type.' )
        df = Io.to_df_watchers(compartments, string, sep=sep)
        return df
    # Adds the infects for infecting compartments.
    def to_df_watchers(
            compartments, # Are infections for valid compartments ? 
            string, # contents, e.g., of input CSV filename for dataframe 
            sep=','): # separator for CSV dataframe
        # Checks the contents of a CSV file.
        try:
            ifh = StringIO(string)
            df = pd.read_csv(ifh, sep=sep,
                    dtype={'passive_watcher':str,'is_difference':int,'compartments':str,'eta_threshold':int})
            if not len(df):
                return None
        except:
            raise ValueError( f'Error: Some data has an incorrect type.' )
        Io.check_watchers(compartments, df)
        return df
    # Checks the watchers for the epidemic.
    def check_watchers(compartments, df): # df describes the active watcher for ending the epidemic.
        passive_watcher_names = [w['name'] for w in to_passive_watchers()] 
        for name in list(df['passive_watcher']):
            if name not in passive_watcher_names:
                raise ValueError(f'Error : "name = "{name}" is not in {passive_watcher_names}.')
        # ? Is the active watcher looking at new cases (is_difference=True) or total cases (is_difference=False) ?
        for i in range(len(df)): 
            cols = ('is_difference','compartments','eta_threshold')
            if all([df.loc[i,col] == '' for col in cols]):
                continue
            is_difference = int(df.loc[i,'is_difference'])
            if not isinstance(is_difference, int) or is_difference not in (0,1):
                raise ValueError(f'is_difference = {is_difference} is not either 0 or 1')
            # compartments separated by ':' for aummation (e.g., I and R) possibly before differencing to discover threshold
            compartments_c = df.loc[i,'compartments']
            compartments0 = compartments_c.split(sep='+')
            if len(compartments0) != len(set(compartments0)):
                raise ValueError(f'Each compartment in {compartments} must be unique.')
            if len(set(compartments0)-set([c for c in compartments])):
                raise ValueError(f'Each compartment in {compartments} must be a valid compartment.')
            # eta_threshold for the active watcher
            eta = int(df.loc[i,'eta_threshold'])
            if not isinstance(eta, int) or not 0 <= eta:
                raise ValueError(f'eta_threshold = {eta} is not an integer 0 <= {eta} = eta')
    # Demonstrates the sequence of CSV files to load the DataFrame-s.        
    def test_df(): 
        DIR = 'Test_Data/'
        #
        dfs_initialize_from = Io.dfs_initialize_from(DIR+'0_initialize_from_0.csv')[0]
        compartments = list(dfs_initialize_from.columns)
        df_add_decay_to = Io.df_add_decay_to(compartments, DIR+'1_add_decay_to_0.csv')
        decay_froms = list(df_add_decay_to['from'])
        df_add_infect_to = Io.df_add_infect_to(compartments, decay_froms, DIR+'2_add_infect_to_0.csv')
        if df_add_infect_to is None:
            pass
        df = Io.df_watchers(compartments, DIR+'3_to_epidemic_with_watchers_0.csv')
        if df is None: # Do nothing with epidemic.
            pass
    #
    # 1 The following routines initialize an epidemic with VALID DataFrames.
    #
    # Initializes an Epidemic from a dataframe. 
    #   Compartment names must contain only letters & underscores.
    #   Assumes check_initialize_from(df) has been called.
    #   i=0 is the 1st column of counts after the column of compartments.
    def initialize_from(df, 
            DELTA=None): # (continuous-time) float for time-grid of the epidemic simulation.
        compartment_names = list(df.columns)
        counts = df.iloc[0,:].values.tolist()
        c2n = dict(zip(compartment_names, counts))
        compartment0counts = [(k,v) for k,v in c2n.items()]
        Compartment.initialize_from(DELTA, *compartment0counts)
    
    def test_initialize_from(): 
        compartment_names = ['s', 'e', 'i', 'r', 'i0', 'r0']
        countss = [[100,1,0,0,0,0]]
        df = pd.DataFrame(countss, columns=compartment_names)
        Io.initialize_from(df)
    # Adds the decays for unstable compartments.
    #   Assumes check_add_decay_to(df) has been called.
    def add_decay_to(df): # dataframe with columns ['from','to','mean','dispersion','prob']
        for i, row in df.iterrows():
            if pd.isna(df.at[i,'dispersion']):
                dispersion = None
            else:
                dispersion = df.at[i,'dispersion']
            Compartment.add_decay_to(
                    name_from=df.at[i,'from'],name_to=df.at[i,'to'], 
                    mean_recovery_time=df.at[i,'mean'], dispersion=dispersion, 
                    prob=df.at[i,'prob'])
    
    def test_add_decay_to(): 
        cols = ['from','to','mean','dispersion','prob']
        lst = [['e','i',2.5,4.0,0.25],['e','i0',2.5,4.0,0.75],['i','r',5.5,0.3,1.0],['i0','r0',5.5,0.3,1.0]]
        df = pd.DataFrame(lst, columns=cols)
        Io.add_decay_to(df)
    # Adds the decays for unstable compartments.
    #   Assumes check_add_infect_to(df) has been called.
    def add_infect_to(
            df, # dataframe with columns ['from','to','infecting','R_0']
            NU=1.0): # multiplier for non-pharmaceutic interventions like school closure
        # The contents of the dataframe are valid.
        for i, row in df.iterrows():
            infecting = df.at[i,'infecting']
            index = Compartment.name2index[infecting]
            compartment = Compartment.compartments[index]
            mean_recovery_time = compartment.decays[0].mean_recovery_time # trivial mixture
            beta = df.at[i,'R_0']/mean_recovery_time
            Compartment.add_infect_to(name_infectious=infecting, beta=beta, name_from=df.at[i,'from'], name_to=df.at[i,'to'])
        Compartment.wrapup(NU)
        
    def test_add_infect_to(): 
        cols = ['from','to','infecting','R_0']
        lst = [['s','e','i',2.0],['s','e','i0',1.5]] # i0 is 0.75 as infectious as i
        df = pd.DataFrame(lst, columns=cols)
        Io.add_infect_to(df, NU=1.0)
    # Initializes the passive and active watchers.
    def to_watchers(df=None): # df describes the active watchers for the epidemic.
        watchers = to_passive_watchers()
        if df is not None:
            # Adds the appropriate active watcher.
            for i in range(len(df)):
                # the relevant passive watcher
                name = df.loc[i,'passive_watcher']
                passive_watcher = next((w for w in watchers if w['name'] == name), None)
                if passive_watcher == None:
                    raise ValueError('No passive watcher is named "{name}".')
                is_difference = int(df.loc[i,'is_difference'])
                eta = int(df.loc[i,'eta_threshold'])
                compartments_c = df.loc[i,'compartments']
                compartments = compartments_c.split(sep='+')
                if is_difference:
                    active_watcher = active_for_difference(passive_watcher, eta, compartments)
                else:
                    active_watcher = active_for_total(passive_watcher, eta, compartments)
                watchers.append(active_watcher)
        return watchers
    # Initializes the epidemic with the watchers.
    #   Adds an active watcher if desired.
    #   Assumes that initialize_from(), add_decay_to(), add_infect_to() have been called.
    #   Deep copy initial epidemic and watchers before epidemic.evolve() to avoid initialization costs.
    def to_epidemic_with_watchers(watchers): 
        epidemic = Epidemic(Compartment.compartments) 
        callback = Callback(watchers)
        epidemic.initialize(callback)
        return epidemic
    #
    # 1 The following routine initialize an epidemic with VALID DataFrames.
    #   It assumes the DataFrames are validated with the check... routines above.
    #
    def test_to_epidemic_with_watchers():         
        cols = ['passive_watcher','is_difference','compartments','eta_threshold']
        lsts = [
            ['all_times',0,'i+r',29],
            ['integer_times',0,'i+r',39],
            ['integer_times',1,'i+r',49]
        ]
        df = pd.DataFrame(lsts, columns=cols) # The df has a single row.
        watchers = Io.to_watchers(df)
        epidemic = Io.to_epidemic_with_watchers(watchers)
        if epidemic is None: # Do nothing with epidemic.
            pass
    #
    # 0 The following routine initialize an epidemic with VALID DataFrames.
    #
    # Demonstrates the sequence of calls required to initialize the epidemic with DataFrame-s.        
    def test(): 
        # The columns of initialize_from give the compartment names.
        compartment_names = ['s', 'e', 'i', 'r', 'i0', 'r0']
        countss = [[100,1,0,0,0,0], [100,2,0,0,0,0], [100,3,0,0,0,0]]
        dfs_initialize_from = pd.DataFrame(countss, columns=compartment_names)
        # 
        cols = ['from','to','mean','dispersion','prob']
        lst = [['e','i',2.5,4.0,0.25],['e','i0',2.5,4.0,0.75],['i','r',5.5,0.3,1.0],['i0','r0',5.5,0.3,1.0]]
        df_add_decay_to = pd.DataFrame(lst, columns=cols)
        #
        cols = ['from','to','infecting','R_0']
        lst = [['s','e','i',2.0],['s','e','i0',1.5]] # i0 is 0.75 as infectious as i
        df_add_infect_to = pd.DataFrame(lst, columns=cols)
        #
        cols = ['passive_watcher','is_difference','compartments','eta_threshold']
        lsts = [
            ['all_times',0,'i+r',29],
            ['integer_times',0,'i+r',39],
            ['integer_times',1,'i+r',49]
        ]
        df_watchers = pd.DataFrame(lsts, columns=cols) # The df has a single row
        # Sets up the Epidemic implicitly with VALID DataFrame-s.
        Io.initialize_from(dfs_initialize_from, DELTA=None)
        Io.add_decay_to(df_add_decay_to)
        Io.add_infect_to(df_add_infect_to, NU=1.0)
        watchers = Io.to_watchers(df_watchers)
        epidemic = Io.to_epidemic_with_watchers(watchers)
        if epidemic is None: # Do nothing with epidemic.
            pass

def main():
    Io.test_df() # Shows how to input and check DataFrame-s.
    Io.test() # Shows how to set up an Epidemic with VALID DataFrame-s.
    
if __name__ == "__main__":
    main()
