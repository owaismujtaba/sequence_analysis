import os

CUR_DIR = os.getcwd() + '/'
CSV_FILES = CUR_DIR +'CSV_WITHOUT_DUPLICATES/'
SPLIT_CSV = CUR_DIR + 'SPLIT_CSV/'


ALLIGN_MATCHED = CUR_DIR + 'Results/ALLIGN_MATCHED/'
COMBINED_ALLIGNED = CUR_DIR + 'Results/COMBINED_ALLIGNED/'
DUPLICATE_INFO = CUR_DIR + 'Results/DUP_INFO/'
SPLIT_ALLIGNED = CUR_DIR + 'Results/SPLIT_ALLIGNED/'
TPRIME_END = CUR_DIR + 'Results/3PRIME_END/'
TPRIME_ALLIGNED = CUR_DIR + 'Results/3PRIME_ALLIGNED/'
EXCEL = TPRIME_ALLIGNED + 'Excel/'

PATTERN = 'Tggctcctcatagggggctcgaacccctgacatttcggttaaaagccgaacgctctaccaactgagctaacaagga'
PATTERN = PATTERN[::-1]
PATTERN = PATTERN.upper()
print('Original Pattern: ', PATTERN)
#PATTERN = PATTERN[30:]