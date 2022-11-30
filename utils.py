import config
import glob
import pdb
import pandas as pd
import os

def get_all_files_in_directory():
    
    
    directory = config.SPLIT_ALLIGNED
    sub_dirs = [os.path.join(directory, x) for x in os.listdir(directory)]
    #pdb.set_trace()
    for folder in sub_dirs:
        if '.ipy' in folder:
            continue
        
        files = next(os.walk(folder), (None, None, []))[2]
        print('Processing {} folder'.format(folder))
        for index in range(len(files)):
            print('Processing {} file'.format(files[index]))
            filename = os.path.join(folder, files[index])
            if index == 0:
                
                full_data = pd.read_csv(filename)
    
            else:
                temp_data = pd.read_csv(filename)
                full_data = pd.concat([full_data, temp_data])
        
        full_data.reset_index(drop=True)
        full_data.drop('Unnamed: 0', axis=1, inplace=True)
        #pdb.set_trace()
        filename = os.path.join(config.COMBINED_ALLIGNED, folder.split('/')[-1]) + '.csv'
        
        full_data.to_csv(filename)
        
    
    
    
    
def clean_combined_matched_files():
    
    allowed_mismatches = 9
    ref_begin = 30
    mismatch_level = 'ML' +str(allowed_mismatches)+'/ML'+str(allowed_mismatches)+'_'
    #pdb.set_trace()
    
    directory = config.COMBINED_ALLIGNED
    files = next(os.walk(directory), (None, None, []))[2]
    index = 0
    for file in files:
        print("Processing {} file".format(file))       
        filename = config.COMBINED_ALLIGNED + file
        data = pd.read_csv(filename)
        print("No of Records: ", data.shape[0])
        seq_length = len(data.sequence[0])
        
        match_threshold = seq_length - allowed_mismatches
        
        
        #pdb.set_trace()
        
        data = data[data.score>=match_threshold]
        data = data[data.ref_start >= ref_begin]
        data.drop('Unnamed: 0', axis=1, inplace=True)
        data.reset_index(drop=True)
        filename = config.TPRIME_END + mismatch_level +file
        #pdb.set_trace()
        data.to_csv(filename)
        print('{} sequences found'.format(data.shape[0]))
        if index == 0:
            full_data = data
        else:
            full_data = pd.concat([full_data, data])
        index += 1
    filename = config.TPRIME_END +mismatch_level + 'full_data.csv'
    print("A total of {} sequences".format(full_data.shape[0]))
    full_data.reset_index(drop=True)
    
    full_data.to_csv(filename) 
    #pdb.set_trace()
    #data = data[data.match_length>=match_threshold]
    
    
    
    
    
clean_combined_matched_files()
#get_all_files_in_directory()   
    