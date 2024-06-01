import os, sys
from probe_design_adapted import designHCR3Probes
from utils import select_output_directory, get_six_digit_date_today


############################# CHANGE SETTINGS BEFORE EACH RUN ##############################
# Please refer to README.md for detailed explanations
GENE_NAMES = [

    ]
GENE_IDS = [

]
GENE_SEQS = [
    
]
HAIRPIN_IDS = [

]
PRB_LENGTH = 20
GC_RANGE = [40, 60]
PRB_SPACING = 2
DG_THRESHOLD = -9
############################################################################################


date = get_six_digit_date_today()
result_path = os.path.join(select_output_directory(), f"HCR3_probe_design_files_{date}")
try:
    os.mkdir(result_path)
except:
    print(f"There is already a directory at {result_path}. Please ensure your selected directory is unique before proceeding.")
    sys.exit(0)

if GENE_SEQS:
    for i in range(len(GENE_NAMES)):
        resultdf = designHCR3Probes(gene_name=GENE_NAMES[i], 
                                    email='czhangyx@berkeley.edu',
                                    sequence=GENE_SEQS[i],
                                    hairpin_id=HAIRPIN_IDS[i], 
                                    db=os.path.join(os.getcwd(), "data/mouse/mouse_refseq_rna"),
                                    result_path=result_path,
                                    prb_length=PRB_LENGTH,
                                    gc_range=GC_RANGE,
                                    prb_space=PRB_SPACING,
                                    dg_thresh=DG_THRESHOLD,
                                    to_excel=True)
else:
    for i in range(len(GENE_NAMES)):
        resultdf = designHCR3Probes(gene_id=GENE_IDS[i], 
                                    gene_name=GENE_NAMES[i], 
                                    email='czhangyx@berkeley.edu',
                                    hairpin_id=HAIRPIN_IDS[i], 
                                    db=os.path.join(os.getcwd(), "data/mouse/mouse_refseq_rna"),
                                    result_path=result_path,
                                    prb_length=PRB_LENGTH,
                                    gc_range=GC_RANGE,
                                    prb_space=PRB_SPACING,
                                    dg_thresh=DG_THRESHOLD,
                                    to_excel=True)