## Dependencies:
Python packages: NumPy, Pandas, Biopython, primer3-py, openpyxl  
Command line tools: Bowtie2  
Other softwares: [Nupack](https://www.nupack.org)

## Installation

## Functions

### Primer design
Primer generation from sequence: run primer_generator.py  
Primer generation from chromosome coordinate: run generate_primers_from_coordinates.py  
Input file format requirement:  
- Excel file

### HCR probe design
Run HCR3_probe_designer.py

### Reverse translation
Input file format requirement: text file with tab-separated values. The file should have 5 columns and 2 rows, including a header row. The second row should include information in the order of IDT username, IDT password, client ID, and client secret.  
Client ID and client secret can be found on IDT's website (see more at https://www.idtdna.com/pages/tools/apidoc):
1. Click your name on the upper right corner of the web page.
2. Go to "My account" -> "API Access"
3. You can request an API key here. After acquiring a key, use your client ID and client secret, as well as your IDT username and password, to access the API with the document format written above.

#### Citations
[HCR probe design](https://doi.org/10.1038/s41587-022-01648-w)  

#### Data Sources
[USeqFISH/HCR databases](https://doi.org/10.22002/b791z-fzd10)  
[CCDS Database](https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi)
[UCSC Genome Browser](https://www.genome.ucsc.edu)
[IDT Codon Optimization](https://www.idtdna.com/restapi/swagger/docs/v1)
