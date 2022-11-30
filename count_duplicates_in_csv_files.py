import pandas as pd
from os import walk
import time
from Bio import SeqIO
import config
import pdb


def count_duplicates():
    file = config.FASTQ_FILES + 'SRR3164724.fastq'
    filename = config.DUPLICATE + file.split('/')[-1] + '.csv'
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
    
    pdb.set_trace()
        
    file = pd.DataFrame(file)
    file = file.drop_duplicates(subset=['sequence'], ignore_index=True)
    file = file.reset_index(drop=True)
    file.to_csv(filename)
    
count_duplicates()