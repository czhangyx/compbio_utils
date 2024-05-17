import numpy as np
import pandas as pd
from utils import select_input_file, report_progress
from utils import import_IDT_information, reverse_translate

################################# CHANGE SETTING BEFORE EACH RUN #################################
IMPORT_FILE_NAME = 'Candidate_Type_I_TM_proteins_sequence_windows.xlsx'
EXPORT_FILE_NAME = 'Candidate_Type_I_TM_proteins_reverse_translated.xlsx'
NAME_COLUMN_NUMBER = 0  # Column number that contains name, zero-indexed
SEQUENCE_COLUMN_NUMBER = 5  # Column number that contains sequence information, zero-indexed
ORGANISM = 'Mus musculus (mouse)'
##################################################################################################


# Import IDT information and data
IDT_file = select_input_file({("text", ".txt"), ("csv", ".csv")})
IDT_info = import_IDT_information(IDT_file)
df = pd.read_excel(IMPORT_FILE_NAME)
array = df.values
titles = df.columns.to_list()
total = len(array)

print('\nReverse translation started...')
data = []
for row in array:
    name = row[NAME_COLUMN_NUMBER]
    seq = row[SEQUENCE_COLUMN_NUMBER].upper()
    assert all([aa in 'QWERTYIPASDFGHKLCVNM' for aa in seq]), \
        f'Sequence for {name} has invalid amino acid'
    data.append({'Name': {name},
                 'Sequence': seq})
dnas = np.array(reverse_translate(data, ORGANISM, *IDT_info))

print('\nGenerating excel file...')
new_df = pd.DataFrame(np.hstack((array, dnas.reshape(-1, 1))))
new_df.columns = titles + [f'{titles[SEQUENCE_COLUMN_NUMBER]} Reverse Translated']
new_df.to_excel(EXPORT_FILE_NAME, index=False)
print('Complete!\n')
