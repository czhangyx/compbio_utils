import os
import pandas as pd
import numpy as np
from utils import get_six_digit_date_today, select_input_file, select_output_directory, report_progress
from utils import get_sequence_from_coordinate, pick_primer, Primer


################################# CHANGE SETTINGS BEFORE EACH RUN ################################
# Please refer to README.md for detailed explanations
CONVERT_FROM_COORDINATES = False
COORDINATE_COLUMN_NUMBER = 1
SEQUENCE_COLUMN_NUMBER = 2
FLANK_SIZE = 50
TARGET_LENGTH = 500
TIGHT_FLANK = True
LEFT_HOMOLOGY_ARM = 'cttttactacctgactagctGCTAGC '
RIGHT_HOMOLOGY_ARM = 'gccctgacttttatgcccagTTAATTAA '
FORBIDDEN = ['GGAGG', 'TAAGGAG', 'TTTTT', 'AAAAA']
LOW_GC = 40
HIGH_GC = 60
LOW_TM = 50
HIGH_TM = 60
LEN_RANGE = range(22, 25)
##################################################################################################


# Import file
import_file_name = select_input_file([("Excel-1", ".xlsx"), ("Excel-2", ".xls"),
                                      ("Excel-3", ".xlsm"), ("Excel-4", ".xlsb"),
                                      ("Excel-5", ".odf"), ("Excel-6", ".ods"),
                                      ("Excel-7", ".odt"), ("CSV", ".csv"), ("TXT", ".txt")])
df = pd.read_excel(import_file_name)
array = df.values
titles = df.columns.to_list()

# Retrieve sequences for genome if necessary
if CONVERT_FROM_COORDINATES:
    print('\n------------ Retrieving sequences from chromosome coordinates ------------')
    coordinates = [peak[COORDINATE_COLUMN_NUMBER] for peak in array]
    total_seqs = []
    current_percentage = 0
    for i in range(len(array)):
        total_seqs.append(get_sequence_from_coordinate(coordinates[i]))
        report_progress(i+1, len(array))
else:
    total_seqs = [peak[SEQUENCE_COLUMN_NUMBER] for peak in array]


# Select forward and reverse primers
print('\n------------ Selecting primers ------------')
print(f'\nSelection criteria:\n{LEN_RANGE[0]} ≤ length ≤ {LEN_RANGE[-1]};\n{LOW_GC}% ≤ GC% ≤ {HIGH_GC}%;\n{LOW_TM}°C ≤ Tm ≤ {HIGH_TM}°C\n')
primers = []
i = 1
for seq in total_seqs:
    primer = Primer()

    for length in LEN_RANGE:
        if primer.fwd and primer.rev:
            break
        fwd_flank = seq[:FLANK_SIZE]
        rev_flank = seq[TARGET_LENGTH+FLANK_SIZE:]
        while len(fwd_flank) >= length and len(rev_flank) >= length:
            if primer.fwd and primer.rev:
                break
            if TIGHT_FLANK:
                pick_primer(primer, 'fwd', fwd_flank[-length:], FORBIDDEN,
                            LOW_GC, HIGH_GC, LOW_TM, HIGH_TM)
                pick_primer(primer, 'rev', rev_flank[:length], FORBIDDEN,
                            LOW_GC, HIGH_GC, LOW_TM, HIGH_TM)
                fwd_flank.pop()
                rev_flank.pop(0)
            else:
                pick_primer(primer, 'fwd', fwd_flank[:length], FORBIDDEN,
                            LOW_GC, HIGH_GC, LOW_TM, HIGH_TM)
                pick_primer(primer, 'rev', rev_flank[-length:], FORBIDDEN,
                            LOW_GC, HIGH_GC, LOW_TM, HIGH_TM)
                fwd_flank.pop(0)
                rev_flank.pop()

    primers.append(primer)
    primer.message(i)
    i += 1

# Generate excel file with primers and their information
print('\nGenerating excel file...')
new_array = [primer.fwd_info() + primer.rev_info() for primer in primers]
for enhancer in new_array:
    enhancer.append(LEFT_HOMOLOGY_ARM + enhancer[0])
    enhancer.append(RIGHT_HOMOLOGY_ARM + enhancer[4])
new_df = pd.DataFrame(np.hstack((array, np.array(new_array))))
new_df.columns = titles + ['fwd primer', 'fwd length', 'fwd GC%', 'fwd Tm',
                           'rev primer', 'rev length', 'rev GC%', 'rev Tm',
                           'full fwd', 'full rev']
date = get_six_digit_date_today()
filename = os.path.splitext(os.path.basename(import_file_name))[0] + f'_with_primers_{date}.xlsx'
export_file_name = select_output_directory(filename)
new_df.to_excel(export_file_name, index=False)
print('Complete!\n')
