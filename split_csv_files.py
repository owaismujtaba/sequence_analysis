import pandas as pd
import config
import pdb
from os import walk

def get_files_from_directory():
    
    
    directory = config.CSV_FILES
    files_in_dir = next(walk(directory), (None, None, []))[2]
    
    files = []
    
    for file in files_in_dir:
        if file.endswith('.csv') and '16.fastq.csv' in file:
            files.append(file)
    
    files_with_path = [directory + '/' + file for file in files]
    #pdb.set_trace()
    return files_with_path



def split_file():
    
    
    n_records_per_file = 500000
    out_dir = config.SPLIT_CSV
    files = get_files_from_directory()
    #pdb.set_trace()
    for file in files:
        data = pd.read_csv(file)
        
        print(file, data.shape)

        for index in range(0, int(data.shape[0]/n_records_per_file)+1):

            filename = out_dir + file.split('/')[-1].split('.')[0] + '_' + str(index) + '.csv'
            print('Making {} file'.format(filename))
            start_index = index * n_records_per_file
            end_index = start_index + n_records_per_file

            temp_data = data[start_index:end_index]

            temp_data.to_csv(filename)
            #pdb.set_trace()
    
split_file()