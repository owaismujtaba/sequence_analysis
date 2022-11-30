import pandas as pd
import os
import config
import pdb

def get_files_from_folder():
    
    directory = config.TPRIME_ALLIGNED + 'ML10/'
    files = next(os.walk(directory), (None, None, []))[2]
    
    files_with_path = [directory + file for file in files]
    
    return files_with_path
    
    
    
def convert_csv_to_excel():
    folder = config.EXCEL + '/ML10.xlsx'
    files = get_files_from_folder()
    #pdb.set_trace()
    writer = pd.ExcelWriter(folder, engine='xlsxwriter')
    filename = files
    for file in files:
        temp_data = pd.read_csv(file)
        sheetname = file.split('/')[-1].split('.')[0]
        temp_data.to_excel(writer, sheet_name = sheetname)
    
    
    writer.save()
    
    
    
convert_csv_to_excel()