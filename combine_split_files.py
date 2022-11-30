import pandas as pd
import config
from os import walk
import pdb

filename = 'SRR31647'

def get_files_from_directory():
    
    
    directory = config.ALLIGNED_CSV
    files_in_dir = next(walk(directory), (None, None, []))[2]
    
    files = []
    
    for file in files_in_dir:
        if file.endswith('.csv') and file.startswith(filename):
            files.append(file)
    #pdb.set_trace()
    return files
    
    


def combine_files():
    
    intrested_columns = ['name', 'seq_type', 'sequence', 'score', 'match_length','seq_start', 'seq_end', 'ref_start', 'ref_end', 'cigar', 'identity']
    
    files = get_files_from_directory()
    global filename
    for i in range(0, len(files)):
        load_file = config.ALLIGNED_CSV + files[i]
        print(load_file)
        if i == 0:
            data = pd.read_csv(load_file)
            data = data[intrested_columns]
            #pdb.set_trace()
        else:
            df = pd.read_csv(load_file)
            df = df[intrested_columns]
            data = pd.concat([data, df])
            data.reset_index(drop=True)
        
    
    filename = config.COMBINED_ALLIGNED + filename + '.csv'
    data = data.drop_duplicates(subset=['seq'], ignore_index=True)
    data = data.reset_index(drop=True)
    data.to_csv(filename)
    
    



combine_files()