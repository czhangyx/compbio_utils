import os, sys
import requests, json
import pandas as pd
import tkinter as tk
from datetime import date
from tkinter import filedialog
from base64 import b64encode
from urllib import request, parse
from Bio.SeqUtils import GC
from Bio.SeqUtils.MeltingTemp import Tm_GC, Tm_NN, Tm_Wallace


"""
TODO
"""
def pick_primer(primer, primer_type, selection, forbidden, low_gc, high_gc, low_tm, high_tm):
    if not any([x in selection for x in forbidden]):
        gc = round(GC(selection), 2)
        tm = round(_Tm(selection), 2)
        if gc in range(low_gc, high_gc+1) and tm in range(low_tm, high_tm+1):
            if not primer.fwd and primer_type == 'fwd':
                primer.update_fwd(selection, gc, tm)
            elif not primer.rev and primer_type == 'rev':
                primer.update_rev(_get_complement(selection), gc, tm)


"""
Prompts user to select an input file.

Input:
    valid_filetypes: a dictionary of tuples of acceptable file type names and extensions.
                     Example: {("text", ".txt"), ("csv", ".csv")}
Returns:
    import_file_name: file path where the program will load data from.
"""
def select_input_file(valid_filetypes: dict) -> str:
    input("""\nYou will be asked to select a file for data input. Press enter to continue.\n""")
    root = tk.Tk()
    root.withdraw()
    import_file_name = filedialog.askopenfilename(title="Please select a file", filetypes=valid_filetypes)
    while not import_file_name:
        os.system("""osascript -e 'tell application "Visual Studio Code" to activate'""")
        input("""\nYou did not select a file. Please press enter to select a file.\n""")
        import_file_name = filedialog.askopenfilename(title="Please select a file", filetypes=valid_filetypes)
    root.destroy()
    return import_file_name


"""
Prompts user to select an output directory.

Input:
    new_dir_name: subdirectory name under the user-selected parent directory.
Returns:
    result_path: result path where the program will generate files into.
"""
def select_output_directory(new_dir_name: str) -> str:
    input("""\nYou will be asked to select a directory for data output. Press enter to continue.\n""")
    root = tk.Tk()
    root.withdraw()
    result_path = filedialog.askdirectory(title="Please select an output directory")
    while not result_path:
        os.system("""osascript -e 'tell application "Visual Studio Code" to activate'""")
        input("""\nYou did not select a directory. Please press enter to select a directory.\n""")
        result_path = filedialog.askdirectory(title="Please select an output directory")
    root.destroy()
    result_path = os.path.join(result_path, new_dir_name)
    try:
        os.mkdir(result_path)
    except:
        print(f"There is already a directory at {result_path}. Please ensure your selected directory is unique before proceeding.")
        sys.exit(0)
    return result_path


"""
Returns a 6-digit representation of today's date.
Throws error if input date_format is not supported.

Input:
    date_format: date representation format, defaulted to mmddyy.
                 Acceptable formats: yymmdd, mmddyy, ddmmyy
"""
def get_six_digit_date_today(date_format: str ='mmddyy') -> str:
    date_format = date_format.lower()
    assert date_format in ['yymmdd', 'mmddyy', 'ddmmyy'], f"Date format of {date_format} is unacceptable"

    year = str(date.today().year % 100)
    month = str(date.today().month)
    day = str(date.today().day)
    if len(month) == 1:
        month = '0' + month
    if len(day) == 1:
        day = '0' + day
    
    if date_format == 'mmddyy':
        return month+day+year
    if date_format == 'yymmdd':
        return year+month+day
    if date_format == 'ddmmyy':
        return day+month+year


"""
Print progress of the current program every milestone.

Inputs:
    count: current iteration count number.
    total: total number of iterations.
    milestone: percentage frequency of reporting.
"""
def report_progress(count: int, total: int, milestone: int=10):
    if (count+1) % (total // milestone) == 0:
        print(f"finished {count+1}/{total}")


"""
Connect to UCSC Genome Browser API to return DNA sequence at specified location.

Inputs:
    coordinate: tuple of coordinate information generated from _process_coordinates.
                Example: (chr3, 1490329, 1491329)
    genome: genome of interest.
"""
def get_sequence_from_coordinate(coordinate: tuple, genome: str ='mm10') -> str:
    assert _genome_in_ucsc(genome), f"Genome {genome} is not supported by data source"

    chromosome, start, end = _process_coordinates(coordinate)
    data = requests.get(f"https://api.genome.ucsc.edu/getData/sequence?genome={genome};chrom={chromosome};start={start};end={end}", timeout=30)
    return data.json()['dna'].upper()


"""
Imports IDT information required for accessing the IDT API.

Input: 
    filepath: the input file path retrieved with select_input_file.
Returns: a list of IDT information in the following order:
    [IDT_username, IDT_password, client_ID, client_secret]
"""
def import_IDT_information(filepath: str) -> list:
    return pd.read_csv(filepath).iloc[0, 0].split('\t')


"""
Reverse translate input peptide sequences to DNA with sequence optimization
using IDT's API.
Please refer to TODO for a list of acceptable organisms.

Inputs:
    data: 
    organism: the target organism for codon optimization.
              Please refer to README.md for a full list of acceptable organisms.
    product_type: the target type of DNA for codon optimization.
                  Acceptable inputs: gblock, gene, megamer
    client_id: IDT client ID.
    client_secret: IDT client secret.
    idt_username: IDT username.
    idt_password: IDT password.
Returns: TODO
"""
def reverse_translate(data, organism: str, product_type: str,
                      IDT_username: str, IDT_password: str, client_ID: str, client_secret: str):
    assert _organism_in_idt(organism), "Chosen organism is not in list"
    assert _product_type_in_idt(product_type), "Chosen product type is not valid"

    access_token = _get_idt_access_token(IDT_username, IDT_password, client_ID, client_secret)
    results = requests.post('https://www.idtdna.com/restapi/v1/CodonOpt/Optimize',
                        json={
                            'organism': organism,
                            'optimizationItems': data,
                            'sequenceType': 'aminoAcid',
                            'productType': product_type
                        },
                        headers={
                            'Authorization': f'Bearer {access_token}',
                            'Content-Type': 'application/json'
                        }).json()
    
    # Clean up result and return a list of DNA sequences
    dnas = []
    for result in results:
        if 'Message' in result:
            dnas.append(None)
        else:
            dnas.append(result['OptResult']['FullSequence'])
    return dnas


class Primer:
    def __init__(self):
        self.fwd = ''
        self.fwd_len = 0
        self.fwd_Tm = None
        self.fwd_gc = None

        self.rev = ''
        self.rev_len = 0
        self.rev_Tm = None
        self.rev_gc = None
    
    def update_fwd(self, fwd, fwd_gc, fwd_Tm):
        self.fwd = fwd
        self.fwd_len = len(fwd)
        self.fwd_Tm = fwd_Tm
        self.fwd_gc = fwd_gc
    
    def update_rev(self, rev, rev_gc, rev_Tm):
        self.rev = rev
        self.rev_len = len(rev)
        self.rev_Tm = rev_Tm
        self.rev_gc = rev_gc
    
    def fwd_info(self):
        return [self.fwd, self.fwd_len, self.fwd_gc, self.fwd_Tm]
    
    def rev_info(self):
        return [self.rev, self.rev_len, self.rev_gc, self.rev_Tm]

    def message(self, i):
        if self.fwd and self.rev:
            print(f'\nSuccessfully picked fwd and rev primers for sequence #{i}')
        elif self.fwd:
            print(f'rev was unsuccessful for sequence #{i}')
        elif self.rev:
            print(f'fwd was unsuccessful for sequence #{i}')
        else:
            print(f'Primer selection was unsuccessful for sequence #{i}')



# HELPER FUNCTIONS
def _Tm(sequence: str) -> float:
    # Calculate Tm from 4 methods and return average
    sequence = sequence.upper()
    assert _is_DNA(sequence), "Input DNA is not valid: contains non AGCT character"

    a, t, g, c = 0, 0, 0, 0
    for base in sequence:
        if base == 'A':
            a += 1
        if base == 'T':
            t += 1
        if base == 'G':
            g += 1
        if base == 'C':
            c += 1
    formula = 64.9 + 41*(g+c-16.4)/(a+t+g+c)
    wallace = Tm_Wallace(sequence)
    gc = Tm_GC(sequence)
    nn = Tm_NN(sequence)
    return (formula+wallace+gc+nn)/4

def _get_complement(sequence: str) -> str:
    # Return complement sequence of input DNA sequence
    sequence = sequence.upper()
    assert _is_DNA(sequence), "Input DNA is not valid: contains non AGCT character"

    result = ''
    for base in sequence:
        if base == 'A':
            result = 'T' + result
        if base == 'T':
            result = 'A' + result
        if base == 'G':
            result = 'C' + result
        if base == 'C':
            result = 'G' + result
    return result

def _is_DNA(sequence: str) -> bool:
    # Return if the input sequence is valid DNA
    return all([base in 'AGCT' for base in sequence])

def _process_coordinates(coordinate: str) -> tuple:
    ## Break coordinates into chromosome, start position, and end position
    chromosome = ""
    chr_done = False
    start = ""
    start_done = False
    end = ""
    while coordinate:
        current = coordinate[0]
        if current == ":":
            chr_done = True
        elif current == "-":
            start_done = True
        elif not chr_done:
            chromosome += current
        elif not start_done:
            start += current
        else:
            end += current
        coordinate = coordinate[1:]
    return (chromosome, int(start)-1, end)

def _genome_in_ucsc(genome: str) -> bool:
    return genome in ['anoCar2', 'bosTau2', 'bosTau3', 'bosTau4', 'bosTau6', 'bosTau7',
                      'bosTau8', 'bosTau9', 'calJac3', 'calJac4', 'canFam1', 'canFam2',
                      'canFam3', 'canFam4', 'canFam5', 'canFam6', 'chlSab2', 'ci3',
                      'danRer10', 'danRer11', 'danRer3', 'danRer4', 'danRer5', 'danRer6',
                      'danRer7', 'equCab1', 'equCab2', 'equCab3', 'fr3', 'galGal2',
                      'galGal3', 'galGal4', 'galGal5', 'galGal6', 'gorGor3', 'gorGor4',
                      'gorGor6', 'hg16', 'hg17', 'hg18', 'hg19', 'hg38', 'hs1', 'macFas5',
                      'melGal1', 'melGal5', 'mm10', 'mm39', 'mm7', 'mm8', 'mm9', 'monDom4',
                      'monDom5', 'nasLar1', 'ornAna1', 'ornAna2', 'oryCun2', 'oryLat2',
                      'oviAri1', 'oviAri3', 'oviAri4', 'panPan2', 'panPan3', 'panTro1',
                      'panTro2', 'panTro3', 'panTro4', 'panTro5', 'panTro6', 'papAnu2',
                      'papAnu4', 'ponAbe2', 'ponAbe3', 'rheMac10', 'rheMac2', 'rheMac3',
                      'rheMac8', 'rn3', 'rn4', 'rn5', 'rn6', 'rn7', 'sacCer1', 'susScr11',
                      'susScr2', 'susScr3', 'taeGut1', 'taeGut2', 'tetNig1', 'tetNig2',
                      'xenTro10', 'xenTro9']

def _get_idt_access_token(IDT_username: str, IDT_password: str,
                          client_ID: str, client_secret: str) -> str:
    authorization_string = b64encode(bytes(client_ID + ":" + client_secret, "utf-8")).decode()
    request_headers = {"Content-Type": "application/x-www-form-urlencoded",
                       "Authorization": "Basic " + authorization_string}
    data_dict = {"grant_type": "password",
                 "scope": "test",
                 "username": IDT_username,
                 "password": IDT_password}
    request_data = parse.urlencode(data_dict).encode()
    post_request = request.Request("https://www.idtdna.com/Identityserver/connect/token", 
                                    data = request_data, 
                                    headers = request_headers,
                                    method = "POST")

    response = request.urlopen(post_request)
    body = response.read().decode()
    
    if (response.status != 200):
        raise RuntimeError("Request failed with error code:" + response.status + "\nBody:\n" + body)
    
    return json.loads(body)["access_token"]

def _organism_in_idt(organism: str) -> bool:
    return organism in ["Drosophila melanogaster",
                        "Escherichia coli K12",
                        "Homo sapiens (human)",
                        "Mus musculus (mouse)",
                        "Pichia pastoris",
                        "Saccharomyces cerevisiae",
                        "Arabidopsis thaliana",
                        "Aspergillus niger",
                        "Azotobacter vinelandii",
                        "Bacillus megaterium",
                        "Bacillus subtilis",
                        "Bifidobacterium longum",
                        "Bombyx mori (silkmoth)",
                        "Bos taurus",
                        "Danio rerio (zebrafish)",
                        "Bradyrhizobium japonicum",
                        "Brassica napus (rape)",
                        "Brevibacillus brevis",
                        "Caenorhabditis elegans (nematode)",
                        "Candida albicans",
                        "Canis familiaris (dog)",
                        "Caulobacter crescentus CB15",
                        "Chlamydia trachomatis D/UW-3/CX",
                        "Chlamydomonas reinhardtii",
                        "Clostridium acetobutylicum ATCC 824",
                        "Corynebacterium glutamicum",
                        "Cricetulus griseus (hamster)",
                        "Cyanophora paradoxa",
                        "Danio rerio",
                        "Dictyostelium discoideum",
                        "Emericella nidulans",
                        "Erwinia carotovora subsp. atroseptica SCRI1043",
                        "Escherichia coli",
                        "Escherichia coli B",
                        "Gallus gallus",
                        "Geobacillus stearothermophilus",
                        "Glycine max (soybean)",
                        "Haemophilus influenzae Rd KW20",
                        "Haloarcula marismortui ATCC 43049 (Halobacterium marismortui)",
                        "Halobacterium salinarum",
                        "Hordeum vulgare subsp vulgare (Barley)",
                        "Klebsiella oxytoca",
                        "Klebsiella pneumoniae",
                        "Kluyveromyces lactis NRRL Y-1140",
                        "Lactobacillus acidophilus",
                        "Lactococcus lactis subsp cremoris",
                        "Leishmania donovani",
                        "Macaca fascicularis",
                        "Magnetospirillum magneticum",
                        "Manduca sexta",
                        "Mannheimia haemolytica",
                        "Medicago sativa",
                        "Methanothermobacter thermautotrophicus str. Delta H",
                        "Moorella thermoacetica",
                        "Mycobacterium tuberculosis H37Rv",
                        "Neisseria gonorrhoeae",
                        "Neurospora crassa",
                        "Nicotiana benthamiana",
                        "Nicotiana tabacum (tobacco)",
                        "Oncorhynchus mykiss (rainbow trout)",
                        "Oryctolagus cuniculus (rabbit)",
                        "Oryza sativa (rice)",
                        "Ovis aries (sheep)",
                        "Petunia x hybrida",
                        "Phaseolus vulgaris (lima bean)",
                        "Pisum sativum (pea)",
                        "Plasmodium falciparum 3D7",
                        "Proteus vulgaris",
                        "Pseudomonas aeruginosa PAO1",
                        "Pseudomonas putida",
                        "Rattus norvegicus (rat)",
                        "Pseudomonas syringae pv tomato str DC3000",
                        "Rhizobium leguminosarum",
                        "Rhodobacter capsulatus",
                        "Rhodobacter sphaeroides",
                        "Salmo salar (Atlantic salmon)",
                        "Salmonella typhimurium LT2",
                        "Schistosoma mansoni",
                        "Schizosaccharomyces pombe",
                        "Schmidtea mediterranea",
                        "Serratia marcescens",
                        "Shewanella frigidimarina",
                        "Simian Virus 40",
                        "Sinorhizobium meliloti 1021",
                        "Solanum lycopersicum (tomato)",
                        "Solanum tuberosum (potato)",
                        "Sorghum bicolor",
                        "Spinacia oleracea (spinach)",
                        "Spodoptera frugiperda",
                        "Staphylococcus aureus subsp. aureus",
                        "Streptococcus mutans UA159",
                        "Streptococcus pneumoniae",
                        "Streptomyces coelicolor",
                        "Strongylocentrotus purpuratus (sea urchin)",
                        "Sus scrofa (pig)",
                        "Synechococcus sp. WH 8102",
                        "Synechococcus sp. PCC 7002",
                        "Synechocystis sp. PCC 6803",
                        "Tetrahymena thermophila",
                        "Thalassiosira pseudonana",
                        "Thermus thermophilus HB8",
                        "Tobacco Mosaic Virus",
                        "Toxoplasma gondii",
                        "Trichoplusia ni",
                        "Triticum aestivum (wheat)",
                        "Trypanosoma brucei",
                        "Trypanosoma cruzi",
                        "Ustilago maydis",
                        "Vaccinia Virus",
                        "Vibrio cholerae O1 biovar eltor str. N16961",
                        "Xenopus laevis",
                        "Yarrowia lipolytica",
                        "Yersinia enterocolitica",
                        "Yersinia pestis",
                        "Zea mays"
                        ]

def _product_type_in_idt(product_type: str) -> bool:
    return product_type in ["gblock", "gene", "megamer"]
