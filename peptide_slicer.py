import pandas as pd
import numpy as np

# CHANGEABLE SECTION
IMPORT_FILE_NAME = 'Candidate_Type_I_TM_proteins_sequence.xlsx'
EXPORT_FILE_NAME = 'Candidate_Type_I_TM_proteins_sequence_windows.xlsx'
MAX_LEN = 50
OVERLAP = 10

# Import data
excel = pd.read_excel(IMPORT_FILE_NAME)
excel_array = np.array(excel)

count = 0
domain_info = []
aa_sequences = []
for domain in excel_array:
    if domain[4] >= 30:
        domain_info.append(domain)
        aa_seq = domain[3]
        front = 0
        seq_windows = ''
        while aa_seq != '':
            seq_windows += (aa_seq[:MAX_LEN]+' ')
            count += 1
            if len(aa_seq) < 50:
                break
            front += OVERLAP
            aa_seq = aa_seq[front:]
        aa_sequences.append(seq_windows)
domain_info = np.array(domain_info, dtype=list)
aa_sequences = np.array(aa_sequences, dtype=list)
print(f'Total oligo count: {count}')

# Generate excel file
print('\nGenerating excel file...')
df = pd.DataFrame(np.hstack((domain_info, aa_sequences.reshape(-1, 1))))
df.columns = ['Gene Names', 'Transmembrane Sequence', 'Transmembrane length', 
              'Cytoplasmic Sequence', 'Cytoplasmic length', 'Cytoplasmic Sequence Windows']
df.to_excel(EXPORT_FILE_NAME, index=False)
print('Complete!\n')
