# CompBio Utils: A collection of Python scripts useful for biological research

## Dependencies:
Python packages: NumPy, Pandas, Biopython, primer3-py, openpyxl  
Command line tools: Bowtie2  
Other softwares: [Nupack](https://www.nupack.org)

## Functions
### **Remember to tune the parameters before running each trial!**  
### 1. Primer design
- Input file format requirement: Excel, .csv, or .txt formats are acceptable; sequences need to contain flanking arms already  
- Run primer_generator.py  
- The program will output an excel file containing the information of generated primers. Full forward and reverse primers with homology arms attached can be found in columns named "full fwd" and "full rev".  

| Parameter name | Description |
|:---------|:---------|
| `CONVERT_FROM_COORDINATES` | `True`: Only chromosome coordinates are available in the input file<br>`False`: Full sequences are available in the input file |
| `COORDINATE_COLUMN_NUMBER` | The column number that contains chromosome coordinates in the input file, zero-indexed (value does not matter if `CONVERT_FROM_COORDINATES` is set to `False`) |
| `SEQUENCE_COLUMN_NUMBER` | The column number that contains sequences in the input file, zero-indexed (value does not matter if `SEQUENCE_COLUMN_NUMBER` is set to `False`) |
| `FLANK_INCLUDED` | `True`: Provided coordinates/sequences include flanking arms for primer search<br>`False`: Provided coordinates/sequences do not include flanking arms for primer search |
| `FLANK_SIZE` | Flanking arm length |
| `TARGET_LENGTH` | Gene length |
| `TIGHT_FLANK` | `True`: Algorithm should remove as much of the flanking arms as possible<br>`False`: Algorithm should retain as much of the flanking arms as possible |
| `LEFT_HOMOLOGY_ARM` | Forward homology arm sequence |
| `RIGHT_HOMOLOGY_ARM` | Reverse homology arm sequence |
| `FORBIDDEN` | A list of motifs primers should not contain |
| `LOW_GC` | Lowest acceptable GC content% |
| `HIGH_GC` | Highest acceptable GC content% |
| `LOW_TM` | Lowest acceptable melting temperature |
| `HIGH_TM` | Highest acceptable melting temperature |
| `LEN_RANGE` | A range of acceptable primer length, excluding homology arms |  


### 2. HCR probe design
### **Remember to tune the parameters before running each trial!**  
- Run HCR3_probe_designer.py
- The program will output an excel file containing the information of generated HCR probes. Full probes can be found in columns named "probe A" and "probe B".

| Parameter name | Description |
|:---------|:---------|
| `GENE_NAMES` | A list of gene names |
| `GENE_IDS` | A list of gene IDs, ordered the same as `GENE_NAMES` (use [NCBI nucleotide index](https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi)) (contents does not matter if using `GENE_SEQS`) |
| `GENE_SEQS` | A list of gene sequences, ordered the same as `GENE_NAMES` (contents does not matter if using `GENE_IDS`) |
| `HAIRPIN_IDS` | A list of insulator IDs, ordered the same as `GENE_NAMES` (refer to [this paper](https://doi.org/10.1038/s41587-022-01648-w)) |
| `PRB_LENGTH` | Desired probe length |
| `GC_RANGE` | A list containing lowest and highest acceptable GC content%, respectively |
| `PRB_SPACING` | Minimum spacing between a probe pair |
| `DG_THRESHOLD` | Lowest acceptable delta G (Gibbs free energy) of a probe (cal/mol) |  


### 3. Reverse translation
*This function requires Internet connection.*
### **Remember to tune the parameters before running each trial!**  
- Input file format requirement:  
    - IDT information: text file with tab-separated values. The file should have 5 columns and 2 rows (including a header row). The second row should include information in the order of IDT username, IDT password, client ID, and client secret. See more on [IDT's website](https://www.idtdna.com/pages/tools/apidoc).  
    - Sequences to be reverse translated: Excel, .csv, or .txt formats are acceptable
- Run reverse_translator.py  
- The program will output an excel file containing the information of reverse translated sequences codon optimized to the selected target organism. 
- Acceptable target organisms (please enter the full name between quotes, including parenthese):
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

| Parameter name | Description |
|:---------|:---------|
| `NAME_COLUMN_NUMBER` | The column number that contains peptide names in the input file, zero-indexed |
| `SEQUENCE_COLUMN_NUMBER` | The column number that contains sequences in the input file, zero-indexed |
| `ORGANISM` | Target organism selected from the list above |


### Please do not modify the location and contents of the following files or the scripts can break!
- utils.py  
- any file in the "data" folder

#### References
[HCR probe design](https://doi.org/10.1038/s41587-022-01648-w)  
[USeqFISH/HCR databases](https://doi.org/10.22002/b791z-fzd10)  
[CCDS Database](https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi)  
[UCSC Genome Browser](https://www.genome.ucsc.edu)  
[IDT Codon Optimization](https://www.idtdna.com/restapi/swagger/docs/v1)  

Please report bugs to czhangyx at berkeley dot edu.
