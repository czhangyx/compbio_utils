import numpy as np
import pandas as pd
import requests
from utils import get_idt_access_token

# CHANGEABLE SECTION
IMPORT_FILE_NAME = 'Candidate_Type_I_TM_proteins_sequence_windows.xlsx'
EXPORT_FILE_NAME = 'Candidate_Type_I_TM_proteins_reverse_translated.xlsx'

# IDT API Information
IDT_username = 'czhangyx'
IDT_password = 'czhangyx037'
client_ID = 'czhangyx'
client_description = 'Yixin Zhang'
client_secret = '38645523-6e00-4897-ab5c-80c97e39e26d'
URL = {
    'organisms': 'https://www.idtdna.com/restapi/v1/CodonOpt/Organisms',
    'sequence types': 'https://www.idtdna.com/restapi/v1/CodonOpt/SequenceTypes',
    'product types': 'https://www.idtdna.com/restapi/v1/CodonOpt/ProductTypes',
    'output': 'https://www.idtdna.com/restapi/v1/CodonOpt/Optimize'
}

# Import data
excel = pd.read_excel(IMPORT_FILE_NAME)
excel_array = np.array(excel)

TOTAL = len(excel_array)

# Helper function
def reverse_translate(data):
    # Connect to IDT codon optimization API (accessible at https://www.idtdna.com/restapi/swagger/docs/v1)
    access_token = get_idt_access_token(client_ID, client_secret, IDT_username, IDT_password)
    organism = requests.get(URL['organisms'], headers={
        'Authorization': f'Bearer {access_token}',
        'Content-Type': 'application/json'
    }).json()[2]
    sequence_type = requests.get(URL['sequence types'], headers={
        'Authorization': f'Bearer {access_token}',
        'Content-Type': 'application/json'
    }).json()[1]
    product_type = requests.get(URL['product types'], headers={
        'Authorization': f'Bearer {access_token}',
        'Content-Type': 'application/json'
    }).json()[0]

    results = requests.post(URL['output'],
                        json={
                            'organism': organism,
                            'optimizationItems': data,
                            'sequenceType': sequence_type,
                            'productType': product_type},
                        headers={
                            'Authorization': f'Bearer {access_token}',
                            'Content-Type': 'application/json'}).json()
    
    # Clean up result and return a string of DNA sequences, separated by comma
    dna = ''
    for result in results:
        window_dna = result['OptResult']['FullSequence']
        print(window_dna)
        dna += (window_dna + ', ')
    return dna

print('\nReverse translation started...')
dnas = []
count = 0
for domain in excel_array:
    data = []
    name = domain[0]
    number = 1
    aa = ''
    for letter in domain[5]:
        if letter == ' ':
            data.append({
                'Name': f'{name} #{number}',
                'Sequence': aa
            })
            aa = ''
            number += 1
        else:
            aa += letter
    dnas.append(reverse_translate(data))
    count += 1
    print(f'{count}/{TOTAL} done\n')
dnas = np.array(dnas)
print(dnas)

# Generate excel file
print('\nGenerating excel file...')
df = pd.DataFrame(np.hstack((excel_array, dnas.reshape(-1, 1))))
df.columns = ['Gene Names', 'Transmembrane Sequence', 'Transmembrane length', 
              'Cytoplasmic Sequence', 'Cytoplasmic length',
              'Cytoplasmic Sequence Windows', 'Cytoplasmic Windows Translated']
df.to_excel(EXPORT_FILE_NAME, index=False)
print('Complete!\n')
