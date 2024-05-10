from probe_design_adapted import designHCR3Probes
import os

# CHANGE SETTING BEFORE EACH RUN
# Obtain gene information from https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi
gene_ids = [ # Use NCBI nucleotide index
    'NM_138942.3', # DBH
]
gene_names = [
    'DBH',
    ]
gene_seqs = [ # Leave as empty array if using gene ids
]
hairpin_ids = [1] # Insulator numbers

prb_length = 24
gc_range = [40, 60]
prb_space = 2
dg_thresh = -9
DATE = '051024'


RESULT_PATH = os.path.join(os.getcwd(), f"HCR3_probe_design_files_{DATE}")
if not os.path.isdir(RESULT_PATH):
    os.mkdir(RESULT_PATH)

if gene_seqs:
    for i in range(len(gene_names)):
        resultdf = designHCR3Probes(#gene_id=gene_ids[i], 
                                    gene_name=gene_names[i], 
                                    email='czhangyx@berkeley.edu',
                                    sequence=gene_seqs[i],
                                    hairpin_id=hairpin_ids[i], 
                                    db=os.path.join(os.getcwd(), "../data/mouse/mouse_refseq_rna"),
                                    result_path=RESULT_PATH,
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
                                    db=os.path.join(os.getcwd(), "../data/mouse/mouse_refseq_rna"),
                                    result_path=RESULT_PATH,
                                    prb_length=prb_length,
                                    gc_range=gc_range,
                                    prb_space=prb_space,
                                    dg_thresh=dg_thresh,
                                    to_excel=True)