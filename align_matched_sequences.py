import pandas as pd
import os
import pdb
import config
import swalign

pattern = config.PATTERN
match = 1
mismatch = -1
scoring = swalign.NucleotideScoringMatrix(match, mismatch)
aligner = swalign.LocalAlignment(scoring)
coding = {'A':'T','T':'A','G':'C','C':'G', 'N':'N'}


def align_all():
    
    files = [x for x in range(17, 18)]
    
    prefix = config.TPRIME_END + '/ML9/ML9_'
    postfix = '.csv'
    
    files = [prefix + str(file) + postfix for file in files]

    
    for file in files:
        
        #pdb.set_trace()
        align_matched_sequences(file)
    
        
    

def align_matched_sequences(file):
    #global pattern
    #file = config.TPRIME_END + '/ML2/ML2_14.csv'
    #file = config.SPLIT_ALLIGNED + 'SRR3164713_0.csv'
    print('Aligning {}'.format(file))
    data = pd.read_csv(file)
    print('Shape', data.shape)
    data.drop('Unnamed: 0', axis=1, inplace=True)
    
    filename = config.TPRIME_ALLIGNED +'/ML9/'+ file.split('/')[-1] 
    #pdb.set_trace()
    pat = config.PATTERN
    pattern = config.PATTERN
    pattern = [char for char in pattern]
    

    initial_append = 10
    dummy = [' ' for i in range(initial_append)]
    dummy_len = len(dummy)
    pattern = dummy + pattern
    
    #pdb.set_trace()
    
    alligned_sequences = []
    alligned_sequences2 = []
    alligned_sequences.append(pattern)
    
    count = 0
    reference_sequenced = []
    for index in data.index:
        new_seq = [' '  for i in range(56)]
        seq_start = int(data.seq_start[index]-1)
        seq_end = int(data.seq_end[index])
        ref_start = int(data.ref_start[index]-1)
        ref_end = int(data.ref_end[index])
        flag = 0
        
        
        
        sequence = data.sequence[index]
        seq_duplication = get_duplicate_count(sequence, file, data.seq_type[index])
        
        print(" Seq {} duplication {}".format(index, seq_duplication))
        
        #pdb.set_trace()
        prefix = sequence[:seq_start]
        postfix = sequence[seq_end:]
        result = aligner.align(pat.replace(' ', ''), sequence)
        qfragment, rfragment = result.dump1()

        prefix = [x for x in prefix]
        postfix = [x for x in postfix]
        qfragment = [x for x in qfragment]
        
        sequence = prefix + qfragment + postfix
        
        
        begin = dummy_len + ref_start - seq_start
        
        pattern2 = pattern
        
        if '-' in rfragment:
            count = count + 1 
            #pdb.set_trace()
            r_prefix = pattern2[:dummy_len+ref_start]
            r_postfix = pattern[dummy_len+ref_start+ref_end:]
            pattern2 = r_prefix + [x for x in rfragment] + r_postfix
            alligned_sequences2.append(pattern2)
            flag = 1
            
            
        #
        new_seq[begin:] = sequence
        
        new_seq[0] = data.name[index]
        new_seq[1] = data.seq_type[index]
        new_seq[2] = data.score[index]
        new_seq[3] = data.match_length[index]
        new_seq[4] = data.identity[index]
        new_seq[5] = seq_start
        new_seq[6] = seq_end
        new_seq[7] = ref_start
        new_seq[8] = ref_end
        new_seq[9] = seq_duplication
        
        if flag == 0:
            alligned_sequences.append(new_seq)
        else:
            alligned_sequences2.append(new_seq)
        
    print('Count :',count)
    
    
    dataframe = pd.DataFrame().from_records(alligned_sequences)
    dataframe.columns = ['col' + str(i) for i in range(dataframe.shape[1])]
    dataframe = put_content_info_across_cols1(dataframe)
    print(dataframe.head())
    dataframe.to_csv(filename)
    #df2 = pd.DataFrame().from_records(alligned_sequences2)
    #df2.columns = ['col' + str(i) for i in range(df2.shape[1])]
    #df2 = put_content_info_across_cols(df2)
    
    
    # putting additional information in the dataframe
    
    
def get_duplicate_count(sequence, file, seq_type):
    global coding
    
    if seq_type == 'seq_c':
        sequence = ''.join([coding[nuc] for nuc in sequence])
    
    filename = 'SRR31647' + file.split('/')[-1].split('_')[1]
    #pdb.set_trace()
    meta_data = pd.read_csv(config.DUPLICATE_INFO+filename)
    temp = meta_data[meta_data['sequence']==sequence]
    try:
        duplicates = temp.n_records.values[0]
    except:
        duplicates = 1
    #pdb.set_trace()
    #print(duplicates)
    
    
    return duplicates



def put_content_info_across_cols1(dataframe):
    
    
    n_records = dataframe.col9.values[1:].sum()
    
    As = []
    Cs = []
    Ts = []
    Gs = []
    
    columns = dataframe.columns[10:]
    for col in columns:
        temp = dataframe[['col9', col]]
        A = 0
        C = 0
        T = 0
        G = 0
        for index in temp.index:
            
            try:
                value =temp.iloc[index]['col9']
                letter = temp.iloc[index][col]
                if letter == None or letter == ' ' or value == ' ':
                    A = A + 0
                    C = C + 0
                    T = T + 0
                    G = G + 0
                else:

                    if letter == 'A':
                        A = A + value 
                    elif letter == 'C':
                        C = C + value 
                    elif letter == 'T':
                        T = T + value 
                    elif letter == 'G':
                        G = G + value 
            except:
                pass
        
        As.append(A/n_records)
        
        Cs.append(C/n_records)
        Ts.append(T/n_records)
        Gs.append(G/n_records)
    dummy = [0 for i in range(10)]
    As = dummy + As
    Cs = dummy + Cs
    Ts = dummy + Ts
    Gs = dummy + Gs

    #pdb.set_trace()
    countACTG = [As, Cs, Ts, Gs]
    cols = dataframe.columns
    countACTG = pd.DataFrame(countACTG, columns=dataframe.columns)
    
    dataframe = pd.concat([dataframe, countACTG])
    
    
    return dataframe
    

  
    
    

        
    
#align_matched_sequences()   
align_all()