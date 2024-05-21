import os
import numpy as np
import pandas as pd
from utils import select_input_file, select_output_directory, get_six_digit_date_today
from utils import import_IDT_information, reverse_translate

################################ CHANGE SETTINGS BEFORE EACH RUN #################################
# Please refer to README.md for detailed explanations
NAME_COLUMN_NUMBER = 0  # Column number that contains name, zero-indexed
SEQUENCE_COLUMN_NUMBER = 5  # Column number that contains sequence information, zero-indexed
ORGANISM = 'Mus musculus (mouse)'  # Refer to README.md for a full list of acceptable organisms
##################################################################################################


# Import IDT information and data
IDT_file = select_input_file({("text", ".txt"), ("csv", ".csv")})
IDT_info = import_IDT_information(IDT_file)
import_file_name = select_input_file([("Excel-1", ".xlsx"), ("Excel-2", ".xls"),
                                      ("Excel-3", ".xlsm"), ("Excel-4", ".xlsb"),
                                      ("Excel-5", ".odf"), ("Excel-6", ".ods"),
                                      ("Excel-7", ".odt"), ("CSV", ".csv"), ("TXT", ".txt")])
df = pd.read_excel(import_file_name)
array = df.values
titles = df.columns.to_list()

print('\nProcessing imported file...')
data = []
for row in array:
    name = row[NAME_COLUMN_NUMBER]
    seq = row[SEQUENCE_COLUMN_NUMBER].upper()
    assert all([aa in 'QWERTYIPASDFGHKLCVNM' for aa in seq]), \
        f'Sequence for {name} has invalid amino acid'
    data.append({'Name': {name},
                 'Sequence': seq})
print('\nReverse translation started, please wait patiently...\n')
dnas = np.array(reverse_translate(data, ORGANISM, *IDT_info))

print('\nGenerating excel file...')
new_df = pd.DataFrame(np.hstack((array, dnas.reshape(-1, 1))))
new_df.columns = titles + [f'{titles[SEQUENCE_COLUMN_NUMBER]} Reverse Translated']
date = get_six_digit_date_today()
filename = os.path.splitext(os.path.basename(import_file_name))[0] + f'_reverse_translated_{date}.xlsx'
export_file_name = select_output_directory(filename)
new_df.to_excel(export_file_name, index=False)
print('Complete!\n')
