import os
from probe_design_adapted import designHCR3Probes
from utils import select_output_directory, get_six_digit_date_today


############################## CHANGE SETTING BEFORE EACH RUN ##############################
# Obtain gene information from https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi
gene_ids = [ # Use NCBI nucleotide index
    'NM_009377.2', # TH
    #'NM_010020.3', # SLC6A3
    #'NM_009891.2', # CHAT
]
gene_names = [
    'TH'#, 'SLC6A3', 'CHAT'
    ]
gene_seqs = [ # Leave as empty array if using gene ids
]
hairpin_ids = [3] # Insulator numbers, see https://doi.org/10.1038/s41587-022-01648-w

prb_length = 20
gc_range = [40, 60]
prb_space = 2
dg_thresh = -9
############################################################################################


date = get_six_digit_date_today()
result_path = select_output_directory(f"HCR3_probe_design_files_{date}")
if gene_seqs:
    for i in range(len(gene_names)):
        resultdf = designHCR3Probes(#gene_id=gene_ids[i], 
                                    gene_name=gene_names[i], 
                                    email='czhangyx@berkeley.edu',
                                    sequence=gene_seqs[i],
                                    hairpin_id=hairpin_ids[i], 
                                    db=os.path.join(os.getcwd(), "data/mouse/mouse_refseq_rna"),
                                    result_path=result_path,
                                    prb_length=prb_length,
                                    gc_range=gc_range,
                                    prb_space=prb_space,
                                    dg_thresh=dg_thresh,
                                    to_excel=True)
else:
    for i in range(len(gene_names)):
        resultdf = designHCR3Probes(gene_id=gene_ids[i], 
                                    gene_name=gene_names[i], 
                                    email='czhangyx@berkeley.edu',
                                    #sequence=gene_seqs[i],
                                    hairpin_id=hairpin_ids[i], 
                                    db=os.path.join(os.getcwd(), "data/mouse/mouse_refseq_rna"),
                                    result_path=result_path,
                                    prb_length=prb_length,
                                    gc_range=gc_range,
                                    prb_space=prb_space,
                                    dg_thresh=dg_thresh,
                                    to_excel=True)