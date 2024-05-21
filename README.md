## Dependencies:
Python packages: NumPy, Pandas, Biopython, primer3-py, openpyxl  
Command line tools: Bowtie2  
Other softwares: [Nupack](https://www.nupack.org)

## Installation

## Functions
### 1. Primer design
Primer generation from sequence: run primer_generator.py  
Primer generation from chromosome coordinate: run generate_primers_from_coordinates.py  
Input file format requirement:  
- Excel file

### 2. HCR probe design
Run HCR3_probe_designer.py

### 3. Reverse translation
*This function requires Internet connection*
Input file format requirement: text file with tab-separated values. The file should have 5 columns and 2 rows, including a header row. The second row should include information in the order of IDT username, IDT password, client ID, and client secret. Client ID and client secret can be found on IDT's website (see more at https://www.idtdna.com/pages/tools/apidoc).
Acceptable target organisms (please enter the full name between quotes, including parenthese):
#### Commonly used
"Drosophila melanogaster", "Escherichia coli K12", "Homo sapiens (human)", "Mus musculus (mouse)", "Pichia pastoris", "Saccharomyces cerevisiae"
#### A
"Arabidopsis thaliana", "Aspergillus niger", "Azotobacter vinelandii"
#### B
"Bacillus megaterium", "Bacillus subtilis", "Bifidobacterium longum", "Bombyx mori (silkmoth)", "Bos taurus", "Bradyrhizobium japonicum", "Brassica napus (rape)", "Brevibacillus brevis"
#### C
"Caenorhabditis elegans (nematode)", "Candida albicans", "Canis familiaris (dog)", "Caulobacter crescentus CB15", "Chlamydia trachomatis D/UW-3/CX", "Chlamydomonas reinhardtii", "Clostridium acetobutylicum ATCC 824", "Corynebacterium glutamicum", "Cricetulus griseus (hamster)", "Cyanophora paradoxa"
#### D
"Danio rerio (zebrafish)", "Danio rerio", "Dictyostelium discoideum"
#### E
"Emericella nidulans", "Erwinia carotovora subsp. atroseptica SCRI1043", "Escherichia coli", "Escherichia coli B"
#### G
"Gallus gallus", "Geobacillus stearothermophilus", "Glycine max (soybean)"
#### H
"Haemophilus influenzae Rd KW20", "Haloarcula marismortui ATCC 43049 (Halobacterium marismortui)", "Halobacterium salinarum", "Hordeum vulgare subsp vulgare (Barley)"
#### K
"Klebsiella oxytoca", "Klebsiella pneumoniae", "Kluyveromyces lactis NRRL Y-1140"
#### L
"Lactobacillus acidophilus", "Lactococcus lactis subsp cremoris", "Leishmania donovani"
#### M
"Macaca fascicularis", "Magnetospirillum magneticum", "Manduca sexta", "Mannheimia haemolytica", "Medicago sativa", "Methanothermobacter thermautotrophicus str. Delta H", "Moorella thermoacetica", "Mycobacterium tuberculosis H37Rv"
#### N
"Neisseria gonorrhoeae", "Neurospora crassa", "Nicotiana benthamiana", "Nicotiana tabacum (tobacco)"
#### O
"Oncorhynchus mykiss (rainbow trout)", "Oryctolagus cuniculus (rabbit)", "Oryza sativa (rice)", "Ovis aries (sheep)"
#### P
"Petunia x hybrida", "Phaseolus vulgaris (lima bean)", "Pisum sativum (pea)", "Plasmodium falciparum 3D7", "Proteus vulgaris", "Pseudomonas aeruginosa PAO1", "Pseudomonas putida", "Pseudomonas syringae pv tomato str DC3000"
#### R
"Rattus norvegicus (rat)", "Rhizobium leguminosarum", "Rhodobacter capsulatus", "Rhodobacter sphaeroides"
#### S
"Salmo salar (Atlantic salmon)", "Salmonella typhimurium LT2", "Schistosoma mansoni", "Schizosaccharomyces pombe", "Schmidtea mediterranea", "Serratia marcescens", "Shewanella frigidimarina", "Simian Virus 40", "Sinorhizobium meliloti 1021", "Solanum lycopersicum (tomato)", "Solanum tuberosum (potato)", "Sorghum bicolor", "Spinacia oleracea (spinach)", "Spodoptera frugiperda", "Staphylococcus aureus subsp. aureus", "Streptococcus mutans UA159", "Streptococcus pneumoniae", "Streptomyces coelicolor", "Strongylocentrotus purpuratus (sea urchin)", "Sus scrofa (pig)", "Synechococcus sp. WH 8102", "Synechococcus sp. PCC 7002", "Synechocystis sp. PCC 6803"
#### T
"Tetrahymena thermophila", "Thalassiosira pseudonana", "Thermus thermophilus HB8", "Tobacco Mosaic Virus", "Toxoplasma gondii", "Trichoplusia ni", "Triticum aestivum (wheat)", "Trypanosoma brucei", "Trypanosoma cruzi"
#### U
 "Ustilago maydis"
#### V
"Vaccinia Virus", "Vibrio cholerae O1 biovar eltor str. N16961",
#### X
"Xenopus laevis"
#### Y
"Yarrowia lipolytica", "Yersinia enterocolitica", "Yersinia pestis"
#### Z
"Zea mays"

### 
#### Citations
[HCR probe design](https://doi.org/10.1038/s41587-022-01648-w)  

#### Data Sources
[USeqFISH/HCR databases](https://doi.org/10.22002/b791z-fzd10)  
[CCDS Database](https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi)
[UCSC Genome Browser](https://www.genome.ucsc.edu)
[IDT Codon Optimization](https://www.idtdna.com/restapi/swagger/docs/v1)
