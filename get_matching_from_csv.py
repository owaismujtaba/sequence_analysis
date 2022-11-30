import swalign
import time
import pandas as pd
import os
from os import walk
from Bio import SeqIO
import pdb
import config
import pandas as pd
import time

def sequence4x(seq):
    coding = {'A':'T','T':'A','G':'C','C':'G', 'N':'N'}
    rev_seq = seq[::-1]
    comp_seq = ''.join([coding[nuc] for nuc in seq])
    rev_comp_seq = comp_seq[::-1]
    
    return seq, rev_seq, comp_seq, rev_comp_seq


def alignment():
    
    match = 1
    mismatch = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    aligner = swalign.LocalAlignment(scoring)
    threshold = 0.90
    match_length = 800
    start = time.time()
    print("Running Allignment ")
    pattern = config.PATTERN
    print("Pattern :",pattern)
    #pdb.set_trace()
    start =time.time()
        
    file = config.SPLIT_CSV + '23/' + 'SRR3164723_10.csv'
    filename = config.SPLIT_ALLIGNED+file.split('/')[-1]
    
    print("Running File ", file)
    file = pd.read_csv(file)    
    
    #pdb.set_trace()
    print("File Shape", file.shape)
    
    
    seq_type = []
    seq_id = []
    result_seq = []
        
    scores = []
    identity = []
    matches = []
    map_cigars = []
    
    map_seq_start = []
    map_seq_end = []
    map_ref_start = []
    map_ref_end = []
   
    
    matched_count = 0
    
    
    for i in range(file.shape[0]):
        
        if i%100 == 0:
            print('Mapped {} Sequences in {}'.format(i, time.time()-start))
            #if i == 500: 
            #    break
        
        seq = file.iloc[i]['sequence']
        name = file.iloc[i]['name']
        
        #pdb.set_trace()
        
        seq, seq_r, seq_c , seq_rc = sequence4x(seq)
        seq_len = len(seq)
        
        
                
        result = aligner.align(pattern, seq_c)
        score = result.identity
        if score > threshold and seq_len-result.matches < match_length:

            print('found seq_c', seq_c)
            print('Sequence : ', seq)
            result_seq.append(seq_c)
            seq_type.append('seq_c')
            seq_id.append(name)
            map_cigars.append(result.cigar)
            identity.append(score)
            scores.append(result.score)
            matches.append(result.matches)
                    
            mapping = result.get_locations()
            map_seq_start.append(mapping[0])
            map_seq_end.append(mapping[1])
            map_ref_start.append(mapping[2])
            map_ref_end.append(mapping[3])
            # pdb.set_trace()
            continue
            
           
                
        result = aligner.align(pattern, seq_r)
        score = result.identity
        if score > threshold and seq_len-result.matches < match_length:

            print('found seq_r', seq_r)
            print('Sequence : ', seq)
            result_seq.append(seq_r)
            seq_type.append('seq_r')
            seq_id.append(name)
            map_cigars.append(result.cigar)
            identity.append(score)
            scores.append(result.score)
            matches.append(result.matches)
                    
            mapping = result.get_locations()
            map_seq_start.append(mapping[0])
            map_seq_end.append(mapping[1])
            map_ref_start.append(mapping[2])
            map_ref_end.append(mapping[3])
            #pdb.set_trace()
            continue
        
        
        result = aligner.align(pattern, seq_rc)
        score = result.identity
        if score > threshold and seq_len-result.matches < match_length:

            print('found seq_rc', seq_rc)
            print('Sequence : ', seq)
            result_seq.append(seq_rc)
            seq_type.append('seq_rc')
            seq_id.append(name)
            
            map_cigars.append(result.cigar)
            identity.append(score)
            scores.append(result.score)
            matches.append(result.matches)
                    
            mapping = result.get_locations()
            map_seq_start.append(mapping[0])
            map_seq_end.append(mapping[1])
            map_ref_start.append(mapping[2])
            map_ref_end.append(mapping[3])
            #pdb.set_trace()
            continue

           
        
                
        result = aligner.align(pattern, seq)
        score = result.identity
        if score > threshold and seq_len-result.matches < match_length:

            print('found seq', seq)
            print('Sequence : ', seq)
            result_seq.append(seq)
            seq_type.append('seq')
            seq_id.append(name)
            map_cigars.append(result.cigar)
            identity.append(score)
            scores.append(result.score)
            matches.append(result.matches)
                    
            mapping = result.get_locations()
            map_seq_start.append(mapping[0])
            map_seq_end.append(mapping[1])
            map_ref_start.append(mapping[2])
            map_ref_end.append(mapping[3])
            #pdb.set_trace()
            continue

           
    result = {'name':seq_id, 
              'seq_type':seq_type, 
              'sequence':result_seq,
              'score':scores,
              'match_length':matches,
              'seq_start':map_seq_start,
              'seq_end':map_seq_end,
              'ref_start':map_ref_start,
              'ref_end':map_ref_end,
              'cigar': map_cigars,
              'identity': identity
             }
    result = pd.DataFrame(result)
    result.to_csv(filename)
    #pdb.set_trace()
        
alignment()