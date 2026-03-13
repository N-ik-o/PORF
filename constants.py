#constants
STOP_CODONS_FW = ("TAA","TAG","TGA")
STOP_CODONS_RV = ("TTA","CTA","TCA")

"""
codon dictionaries :
    - keys: codon strings (e.g., 'ATG')
    - values: single-letter amino acid codes (e.g., 'M')
      Stop codons are represented as '*'.
    
Notes
    Uses standard genetic code. 
    IF ALTERNATE CODONS ARE USED IN HOST, ALTER CODES ACCORDINGLY!
"""

CODON_2_AA_FW = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L',
    'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S',
    'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 
    'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 
    'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 
    'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 
    'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 
    'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 
    'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 
    'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 
    'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 
    'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 
    'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 
    'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 
    'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 
    'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 
    'GGG': 'G'}

CODON_2_AA_RV = {'TTT': 'K', 'TTC': 'E', 'TTA': '*', 
    'TTG': 'Q', 'TCT': 'R', 'TCC': 'G', 'TCA': '*', 
    'TCG': 'R', 'TAT': 'I', 'TAC': 'V', 'TAA': 'L', 
    'TAG': 'L', 'TGT': 'T', 'TGC': 'A', 'TGA': 'S', 
    'TGG': 'P', 'CTT': 'K', 'CTC': 'E', 'CTA': '*', 
    'CTG': 'Q', 'CCT': 'R', 'CCC': 'G', 'CCA': 'W', 
    'CCG': 'R', 'CAT': 'M', 'CAC': 'V', 'CAA': 'L', 
    'CAG': 'L', 'CGT': 'T', 'CGC': 'A', 'CGA': 'S', 
    'CGG': 'P', 'ATT': 'N', 'ATC': 'D', 'ATA': 'Y', 
    'ATG': 'H', 'ACT': 'S', 'ACC': 'G', 'ACA': 'C', 
    'ACG': 'R', 'AAT': 'I', 'AAC': 'V', 'AAA': 'F', 
    'AAG': 'L', 'AGT': 'T', 'AGC': 'A', 'AGA': 'S', 
    'AGG': 'P', 'GTT': 'N', 'GTC': 'D', 'GTA': 'Y', 
    'GTG': 'H', 'GCT': 'S', 'GCC': 'G', 'GCA': 'C', 
    'GCG': 'R', 'GAT': 'I', 'GAC': 'V', 'GAA': 'F', 
    'GAG': 'L', 'GGT': 'T', 'GGC': 'A', 'GGA': 'S', 
    'GGG': 'P'}
