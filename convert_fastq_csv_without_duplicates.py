import swalign
from os import walk
import config
import swalign
from Bio import SeqIO
import pdb
import pandas as pd
import time



def fastq_to_csv():
    file = config.FASTQ_FILES + 'SRR3164724.fastq'
    filename = config.CSV_FILES + file.split('/')[-1] + '.csv'
    #pdb.set_trace()
    print("Processing", file)
    count = 0
    sequences = []
    ids = []
    start = time.time()
        
        
    for record in SeqIO.parse(file, 'fastq'):
            count = count+1
            if count % 100000 == 0:
                print('Processed {} sequnces time: {}'.format(count, time.time()-start))
                #break
    
            sequences.append(str(record.seq))
            #pdb.set_trace()
            ids.append(record.name)
            
    file = {'name':ids, 'sequence': sequences}
    #pdb.set_trace()
        
    file = pd.DataFrame(file)
    file = file.drop_duplicates(subset=['sequence'], ignore_index=True)
    file = file.reset_index(drop=True)
    file.to_csv(filename)
        
                
        
        

fastq_to_csv()
        
