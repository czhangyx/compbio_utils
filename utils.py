import os, sys
import tkinter as tk
from tkinter import filedialog
from Bio.SeqUtils import GC
from Bio.SeqUtils.MeltingTemp import Tm_GC, Tm_NN, Tm_Wallace

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
Prompts user to select an output directory.
"""
def select_output_directory(new_dir_name):
    root = tk.Tk()
    root.withdraw()
    os.system("""osascript -e 'tell application "Visual Studio Code" to activate'""")
    input("""\nYou will be asked to select a directory for data output. Press enter to continue.\n""")
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


# HELPER FUNCTIONS
## AVOID DIRECTLY CALLING THESE OUTSIDE
def _Tm(sequence: str):
    sequence = sequence.upper()
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

def _get_complement(sequence: str):
    sequence = sequence.upper()
    assert __is_DNA(sequence), "Input DNA is not valid: contains non AGCT character"

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

def __is_DNA(sequence: str):
    return all([base in 'AGCT' for base in sequence])


