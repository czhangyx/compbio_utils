import numpy as np
import pandas as pd

# CHANGEABLE SECTION
IMPORT_FILE_NAME = 'Candidate_Type_I_TM_proteins_reverse_translated.xlsx'
REFERENCE_FILE_NAME = 'SmallerList.xlsx'
EXPORT_FILE_NAME = 'Candidate_Type_I_TM_proteins_reverse_translated_shortened.xlsx'

# Import data
excel_array = np.array(pd.read_excel(IMPORT_FILE_NAME))
ref_array = np.array(pd.read_excel(REFERENCE_FILE_NAME))
domains = []

for domain in excel_array:
    for name in domain[0].split(' '):
        if name in ref_array:
            if domain[4] <= 50:
                domain[-1] = domain[-1][:-2]
                domains.append(np.hstack((np.array([name]), domain[1:])))
            else:
                number = 1
                seq_list = domain[-1].split(', ')
                for seq in seq_list[:-1]:
                    x = np.hstack((np.array([f'{name}_{number}']), domain[1:-1]))
                    x = np.hstack((x, seq))
                    domains.append(x)
                    number += 1
domains = np.array(domains)
domains = domains[domains[:, 0].argsort()]

# Modify sequences
SpeI = 'ACTAGT'
AgeI = 'ACCGGT'
FRONT_HA = 'TTATCACCCTTTACTGCAGA'
BACK_HA = 'GCCACGAACTTCTCTCTGTT'
SG = ['TCC', 'TCG', 'TCA', 'TCT', 'AGT', 'AGC', 'GGC', 'GGG', 'GGA', 'GGT']
for domain in domains:
    spacer = ''
    seq = domain[-1]
    ## Construct spacer
    if len(seq) < 138:
        spacer = SpeI
        while len(spacer) + 6 < 150 - len(seq):
            spacer += np.random.choice(SG)
        spacer += AgeI
    ## Add homology arms and spacer to sequences
    domain[-1] = FRONT_HA + seq + spacer + BACK_HA

# Generate excel file
print('\nGenerating excel file...')
df = pd.DataFrame(domains)
df.columns = ['Gene Names', 'Transmembrane Sequence', 'Transmembrane length', 
              'Cytoplasmic Sequence', 'Cytoplasmic length',
              'Cytoplasmic Sequence Windows', 'Cytoplasmic Windows Translated']
df.to_excel(EXPORT_FILE_NAME, index=False)
print('Complete!\n')
