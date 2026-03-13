import constants as const
from warnings import warn


"""
------------------------
READ INPUT FILE FUNCTION
------------------------
"""
""" Reads a FASTA file and processes each sequence.

    This function reads a FASTA-formatted input file, concatenates
    multi-line sequences into a single string, and passes each
    complete sequence to the `scan_sequence` function along with
    the specified output name.

    The function assumes that sequences in the FASTA file are
    separated by header lines beginning with the ">" character.

    Args:
        input_file (str): Path to the input FASTA file.
        output_name (str): Name used by `scan_sequence`
            to label output results.

    Returns:
        None

    Notes:
        - Each sequence is processed immediately after its header
          line is encountered.
        - The final sequence in the file is processed after the
          loop completes.
    """

def read_fasta(input_file,output_name):
    with open (input_file) as infile:
        sequence_parts = []
        for line in infile:
            current_line = line.strip()
            if current_line.startswith(">"):
                sequence = "".join(sequence_parts)
                scan_sequence(sequence,output_name)
                sequence_parts = []
            else: 
                sequence_parts.append(current_line)
        sequence = "".join(sequence_parts)
        scan_sequence(sequence, output_name)


"""
-----------------------
SCAN SEQUENCE FUNCTIONS
-----------------------
active one defined by argparser
"""
"""
This module provides wrapper functions for scanning nucleotide 
sequences and writing results to different output formats.
The active function is selected via argument parsing (argparser)
in the main program.

Functions:
    scan_sequence_both:
        Scans a sequence and writes results to both FASTA (.fa)
        and GFF (.gff) output files.

    scan_sequence_fasta:
        Scans a sequence and writes results to a FASTA (.fa) file only.

    scan_sequence_gff:
        Scans a sequence and writes results to a GFF (.gff) file only.

All functions share the same parameters:

Args:
    sequence (str): The processed nucleotide sequence to analyze
    (i.e., no FASTA headers and no newline characters).

    output_name (str): Base name for the output file(s), without file
        extension. Files are created in the `output_files/` directory.

Notes:
    - Files are opened in append mode ("a"), so existing content
      is preserved.
    - All wrapper functions call `scan_frames`, which performs
      the core frame-scanning logic.
"""
def scan_sequence_both(sequence, output_name):
    fasta_file = f"output_files/{output_name}.fa"
    gff_file = f"output_files/{output_name}.gff"
    with open(fasta_file, "a") as f:
        with open(gff_file, "a") as g:
            scan_frames(sequence, f, g)

def scan_sequence_fasta(sequence, output_name):
    fasta_file = f"output_files/{output_name}.fa"
    with open(fasta_file, "a") as f:
        scan_frames(sequence, f, f)

def scan_sequence_gff(sequence, output_name):
    gff_file = f"output_files/{output_name}.gff"
    with open(gff_file, "a") as g:
        scan_frames(sequence, g, g)


"""
--------------------
SCAN FRAMES FUNCTION
--------------------
"""
"""
Scans all reading frames of a sequence.

Iterates over the three possible reading frames (frameshifts 0, 1, and 2)
for a given nucleotide sequence and delegates ORF scanning to `scan_frame`.

Args:
    sequence (str):
        The preprocessed nucleotide sequence to scan for open reading frames (ORFs).

    f (file object):
        File handle for FASTA output. Passed to downstream scanning and writer functions.

    g (file object):
        File handle for GFF output. Passed to downstream scanning and writer functions.
"""

def scan_frames(sequence, f, g):
    for frameshift in range(3):
        scan_frame(sequence, frameshift, f, g)

"""
--------------------
SCAN FRAME FUNCTIONS
--------------------
active one defined by argparser
"""
"""
These functions evaluate a single reading frame of a nucleotide sequence
and extend open reading frames (ORFs) according to forward and/or
reverse strand logic. The active scan_frame function is chosen via 
argument parsing.

Functions:
    scan_frame_both(sequence, frameshift, f, g)
        Scans both forward and reverse strands for ORFs in the given
        reading frame. Iteratively calls `scan_strain_fw` and
        `scan_strain_rv` until the end of the sequence is reached.

    scan_frame_fw(sequence, frameshift, f, g)
        Scans only the forward strand for ORFs in the given reading frame.
        Iteratively calls `scan_strain_fw` until the end of the sequence
        is reached.

Args:
    sequence (str):
        The nucleotide sequence to scan.

    frameshift (int):
        The reading frame offset (0, 1, or 2) to begin scanning.

    f (file object):
        File handle for FASTA output. Passed to strain scanning and writer functions.

    g (file object):
        File handle for GFF output. Passed to strain scanning and writer functions.

Notes:
    - ORFs are tracked using boolean flags (e.g., `extend_orf_fw`,
      `extend_orf_rv`) to determine whether an ORF is currently
      being extended.
"""
        
def scan_frame_both(sequence, frameshift, f, g):
    extend_orf_fw = False
    extend_orf_rv = False
    end_fw = 3 + frameshift 
    end_rv = 3 + frameshift
    start_fw = 0
    start_rv = 0
    while end_fw <= len(sequence) or end_rv <= len(sequence):        
        if end_fw <= len(sequence):
            extend_orf_fw, start_fw, end_fw = scan_strain_fw(sequence, extend_orf_fw, start_fw, end_fw, f, g)
        if end_rv <= len(sequence):
            extend_orf_rv, start_rv, end_rv = scan_strain_rv(sequence, extend_orf_rv, start_rv, end_rv, f, g)  

def scan_frame_fw(sequence, frameshift, f, g):
    extend_orf_fw = False
    end_fw = 3 + frameshift 
    start_fw = 0
    while end_fw <= len(sequence):         
        extend_orf_fw, start_fw, end_fw = scan_strain_fw(sequence, extend_orf_fw, start_fw, end_fw, f, g)

"""
---------------------
SCAN STRAIN FUNCTIONS
---------------------
"""
"""
Scan a single strain (forward or reverse) for ORF extension.

These functions handle the ORF scanning within a single reading
frame for either the forward or reverse strand. Depending on whether 
an ORF is currently being extended, the function either initiates a 
search for a start codon (forward) or a stop codon (reverse), 
or extends the existing ORF to the next codon.

Functions:
    scan_strain_fw(sequence, extend_orf_fw, start_fw, end_fw, f, g)
        Processes the forward strand of a reading frame.

    scan_strain_rv(sequence, extend_orf_rv, start_rv, end_rv, f, g)
        Processes the reverse strand of a reading frame.

Args:
    sequence (str):
        The nucleotide sequence to scan.

    extend_orf_fw / extend_orf_rv (bool):
        Indicates whether an ORF is currently being extended.

    start_fw / start_rv (int):
        The start position of the current ORF. Updated when a start
        (forward) or stop (reverse) codon is found.

    end_fw / end_rv (int):
        The end position of the current ORF. 
        Updated as the ORF is extended.

    f (file object):
        File handle for FASTA output. Passed to ORF extension functions.

    g (file object):
        File handle for GFF output. Passed to ORF extension functions.

Returns:
    tuple:
        - Updated ORF extension flag (bool)
        - Updated start position (int)
        - Updated end position (int)
"""

def scan_strain_fw(sequence, extend_orf_fw, start_fw, end_fw, f, g):
    if extend_orf_fw == False: 
        extend_orf_fw, start_fw, end_fw = atg_search(sequence, extend_orf_fw, end_fw)    
    else: 
        extend_orf_fw, end_fw = orf_extension_fw(sequence, extend_orf_fw, start_fw, end_fw, f, g)
    return extend_orf_fw, start_fw, end_fw

def scan_strain_rv(sequence, extend_orf_rv, start_rv, end_rv, f, g):
    if extend_orf_rv == False: 
        extend_orf_rv, start_rv, end_rv = stop_search(sequence, extend_orf_rv, end_rv)
    else: 
        extend_orf_rv, start_rv, end_rv = orf_extension_rv(sequence, extend_orf_rv, start_rv, end_rv, f, g)
    return extend_orf_rv, start_rv, end_rv

"""
----------------------------------
STARTING POSITION SEARCH FUNCTIONS
----------------------------------
"""
"""
Search for ORF start or stop positions in a nucleotide sequence.

These functions locate the starting codon for forward-strand ORFs
(`atg_search`) or the stop codon for reverse-strand ORFs (`stop_search`).
They are used to initiate ORF extension when no ORF is currently active.

Functions:
    atg_search(sequence, extend_orf_fw, end_fw)
        Searches for a forward-strand start codon predifined by Argument Parsing at the 
        current codon position. Updates the ORF extension flag and start/end positions.

    stop_search(sequence, extend_orf_rv, end_rv)
        Searches for a reverse-strand stop codon at the current codon
        position. Updates the ORF extension flag and start/end positions.

Args:
    sequence (str):
        The nucleotide sequence to scan.

    extend_orf_fw / extend_orf_rv (bool):
        Indicates whether an ORF is currently being extended. If False,
        the function performs a start/stop codon search.

    end_fw / end_rv (int):
        The end position of the current codon being checked. Used to
        determine the triplet to evaluate.

Returns:
    tuple:
        - Updated ORF extension flag (bool)
        - Start position of the ORF (int, 0 if not found)
        - Updated end position of the codon (int)
"""


def atg_search(sequence, extend_orf_fw, end_fw):
    if sequence[end_fw - 3 : end_fw] in start_codons: 
        extend_orf_fw = True
        start_fw = end_fw - 3  
    else:
        end_fw += 3
        start_fw = 0
    return extend_orf_fw, start_fw, end_fw

def stop_search(sequence, extend_orf_rv, end_rv):
    if sequence[end_rv - 3 : end_rv] in const.STOP_CODONS_RV:
        extend_orf_rv = True
        start_rv = end_rv - 3 
    else:
        end_rv += 3
        start_rv = 0
    return extend_orf_rv, start_rv, end_rv

"""
-----------------------
ORF EXTENSION FUNCTIONS
-----------------------
"""
"""
Extend open reading frames (ORFs) until a stop codon (forward) or
a reverse-complemented start codon (reverse) is reached. 
The reverse-complemented start codon are defined by Argument Parsing.

These functions continue the extension of an ORF codon-by-codon for
either the forward or reverse strand. After successful completion,
the ORF is finalized and written to the specified output files.

Functions:
    orf_extension_fw(sequence, extend_orf_fw, start_fw, end_fw, f, g)
        Extends an ORF on the forward strand and writes it upon reaching
        a stop codon.

    orf_extension_rv(sequence, extend_orf_rv, start_rv, end_rv, f, g)
        Extends an ORF on the reverse strand and writes it upon reaching
        a start codon, if no new stop codon is reached. If a new stop codon is 
        reached, start of the ORF is set to tis position.

Args:
    sequence (str):
        The nucleotide sequence being scanned.

    extend_orf_fw / extend_orf_rv (bool):
        Indicates whether the ORF is currently being extended. 
        Set to False when the ORF is completed.

    start_fw / start_rv (int):
        Start position of the current ORF. Used for writing the ORF
        once completed.

    end_fw / end_rv (int):
        End position of the current codon. Updated as the ORF extends.

    f (file object):
        File handle for FASTA output. Passed to the writer function.

    g (file object):
        File handle for GFF output. Passed to the writer function.

Returns:
    tuple:
        - Updated ORF extension flag (bool)
        - Updated end position (int)
        - For reverse strand only: updated start position (int)
"""

#forward
def orf_extension_fw(sequence, extend_orf_fw, start_fw, end_fw, f, g): 
    if sequence[end_fw - 3 : end_fw] not in const.STOP_CODONS_FW:
        end_fw += 3
    else: 
        writer(sequence,start_fw,end_fw,f,g,"forward")
        extend_orf_fw = False
    return extend_orf_fw, end_fw 

#reverse
def orf_extension_rv(sequence, extend_orf_rv, start_rv, end_rv, f, g):
    if sequence[end_rv - 3 : end_rv] not in const.STOP_CODONS_RV and sequence[end_rv - 3 : end_rv] not in start_codons_rv:
        end_rv += 3
    elif sequence[end_rv - 3 : end_rv] in const.STOP_CODONS_RV:
        start_rv = end_rv - 3
        end_rv += 3 
    else: 
        writer(sequence,start_rv,end_rv,f,g,"reverse")
        extend_orf_rv= False
    return extend_orf_rv, start_rv, end_rv


"""
-----------------------
OUTPUT WRITER FUNCTIONS
-----------------------
active one defined by argparser
"""
"""
Write predicted ORFs to FASTA and/or GFF output files.

These functions handle the formatting and writing of open reading frames
(ORFs) once they have been identified. Additionally, they verify whether
the ORF length satisfies the minimum cutoff threshold.
Output can be in nucleotide or amino acid FASTA format, GFF annotation 
format, or both simultaneously.
The active writer function is chosen via argument parsing.

Functions:
    writer_fasta_aa(sequence, start_fw, end_fw, f, g, orientation)
        Writes the ORF as an amino acid sequence to a FASTA file.

    writer_fasta_nuc(sequence, start_fw, end_fw, f, g, orientation)
        Writes the ORF as a nucleotide sequence to a FASTA file.

    writer_gff(sequence, start_fw, end_fw, f, g, orientation)
        Writes the ORF location and frame information to a GFF file.

    writer_both_aa(sequence, start_fw, end_fw, f, g, orientation)
        Writes the ORF as an amino acid FASTA sequence and as a GFF entry.

    writer_both_nuc(sequence, start_fw, end_fw, f, g, orientation)
        Writes the ORF as a nucleotide FASTA sequence and as a GFF entry.

Args:
    sequence (str):
        The nucleotide sequence of the ORF to write.

    start_fw (int):
        Start position of the ORF in the sequence.

    end_fw (int):
        End position of the ORF in the sequence.

    f (file object):
        File handle for FASTA output.

    g (file object):
        File handle for GFF output.

    orientation (str):
        Orientation of the ORF; either 'forward' for forward strand or 'reverse'
        for reverse strand.

Global Variables:
    cutoff_value (int):
        Minimum ORF length threshold in nucleotides. ORFs shorter
        than this value are ignored.
"""
def writer_fasta_aa(sequence,start_fw,end_fw,f,g,orientation):
    orf = sequence[start_fw: end_fw]
    if len(orf) >= cutoff_value:
        if orientation == "forward":
            aa_seq = transcribe(codon_2_aa_fw, orf)
        else:
            aa_seq = transcribe(codon_2_aa_rv, orf)
            aa_seq = aa_seq[::-1]
        fasta_header_fw = f">Start: {start_fw}; end: {end_fw}; orientation: {orientation}"
        f.write(fasta_header_fw + "\n" + aa_seq + "\n")

def writer_fasta_nuc(sequence,start_fw,end_fw,f,g,orientation):
    orf = sequence[start_fw: end_fw]
    if len(orf) >= cutoff_value:
        fasta_header_fw = f">Start: {start_fw}; end: {end_fw}; orientation: {orientation}"
        f.write(fasta_header_fw + "\n" + orf + "\n")

def writer_gff(sequence,start_fw,end_fw,f,g,orientation):
    if end_fw - start_fw >= cutoff_value:
        gff_entry = f"SeqLoc\tsource\torf\t{str(start_fw)}\t{str(end_fw)}\tscore\t{orientation}\t{str(start_fw%3)}"
        g.write(gff_entry + "\n")

def writer_both_aa(sequence,start_fw,end_fw,f,g,orientation):
    orf = sequence[start_fw: end_fw]
    if len(orf) >= cutoff_value:
        if orientation == "forward":
            aa_seq = transcribe(codon_2_aa_fw, orf)
        else:
            aa_seq = transcribe(codon_2_aa_rv, orf)
            aa_seq = aa_seq[::-1]
        fasta_header_fw = f">Start: {start_fw}; end: {end_fw}; orientation: {orientation}"
        f.write(fasta_header_fw + "\n" + aa_seq + "\n")
        gff_entry = f"SeqLoc\tsource\torf\t{str(start_fw)}\t{str(end_fw)}\tscore\t{orientation}\t{str(start_fw%3)}"
        g.write(gff_entry + "\n")

def writer_both_nuc(sequence,start_fw,end_fw,f,g,orientation):
    orf = sequence[start_fw: end_fw]
    if len(orf) >= cutoff_value:
        fasta_header_fw = f">Start: {start_fw}; end: {end_fw}; orientation: {orientation}"
        f.write(fasta_header_fw + "\n" + orf + "\n")
        gff_entry = f"SeqLoc\tsource\torf\t{str(start_fw)}\t{str(end_fw)}\tscore\t{orientation}\t{str(start_fw%3)}"
        g.write(gff_entry + "\n")


"""
-------------------
TRANSCRIBE FUNCTION
-------------------
"""
"""
Translates a nucleotide sequence into an amino acid sequence
in the direction specified by the dictionary used as parameter.

Parameters
----------
codon_dictionary : dict
    Dictionary mapping codons (str) to amino acids (str).

sequence : str
    DNA sequence consisting of nucleotide characters (e.g., "ATGCGT...").

Returns
-------
str
    Amino acid sequence translated from the input DNA sequence.
"""

def transcribe(codon_dictionary,sequence):
    i = 0
    aa_sequence = ""
    if len(sequence) % 3 != 0:
        warn("Length not a multiple of 3. Translating until the highest multiple of 3 < len(sequence)")
    while i + 3 <= len(sequence):
        codon = sequence[i:i+3]
        if codon not in codon_dictionary.keys():
            warn("Ambiguous base detected.")
        else: aa_sequence += codon_dictionary[codon]
        i += 3
    return aa_sequence



"""
----------
ARGPARSER
----------
"""
"""
Command-line argument parser for ORF prediction.

This parser defines the command-line interface for running the ORF
prediction pipeline. Users can specify input/output files, output
formats, sequence types, scanning directions, and additional parameters 
(start codons, cutoff value) affecting ORF detection.

Arguments:
    -f, --format (str, optional):
        Output format for ORF predictions. Choices are:
        - "fasta" : write FASTA file only
        - "gff"   : write GFF file only (default)
        - "both"  : write both FASTA and GFF files

    output (str):
        Base name for the output file(s). Do not include file extensions;
        extensions will be added automatically based on the chosen format.

    input (str):
        Path to the input FASTA file containing nucleotide sequences.

    -fo, --fasta_output (str, optional):
        If FASTA output is chosen, specify the sequence type:
        - "aa"  : translate ORFs to amino acid sequences
        - "nuc" : output nucleotide sequences (default)

    -d, --direction (str, optional):
        Scanning direction for ORF detection:
        - "fw"   : forward strand only (default)
        - "both" : scan both forward and reverse strands

    -c, --cutoff (int, optional):
        Minimum ORF length threshold specified in amino acids.
        ORFs shorter than this cutoff value are ignored.
        Default: 1

    -s, --start_codons (str, optional):
        Defines which start codons are considered during ORF detection.
        Choices:
            "ATG"           : canonical start codon only (default)
            "+GTG"          : include GTG in addition to ATG
            "+TTG"          : include TTG in addition to ATG
            "+GTG+TTG"      : include both GTG and TTG in addition to ATG

"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--format", choices = ["fasta", "gff", "both"], default = "gff",
                    help = "Choose the output format(s).")
parser.add_argument("output", help = "Enter your output filename without extension.")
parser.add_argument("input", help =  "Input path to your input fasta file.")
parser.add_argument("-fo", "--fasta_output", choices = ["aa","nuc"], default = "nuc", 
                    help= "Choose fasta output style: amino acids (aa) or nucleotides (nuc).")
parser.add_argument("-d", "--direction", choices=["fw","both"], default = "fw", 
                    help = "Choose scanning directions: forward (fw) or forward and reverse (both)")
parser.add_argument("-c", "--cutoff", default = "1", help = "Choose cutoff value: minimum amino acid length.")
parser.add_argument("-s", "--start_codons", default = "ATG", choices=["ATG","+GTG","+TTG","+GTG+TTG"], 
                     help = "Choose start codons. (ATG), (+GTG), (+TTG), (+GTG+TTG)")

args = parser.parse_args()

"""
------------------------------
ARGPARSE CONDITIONAL FUNCTIONS
------------------------------
"""
"""
Select active ORF prediction functions based on user arguments.

This section maps the command-line arguments to the appropriate
functions for ORF prediction, frame scanning, and output writing.

Behavior:
    - `args.format` determines which prediction function to use:
        - "both"  : predict_orfs_both
        - "fasta" : predict_orfs_fasta
        - "gff"   : predict_orfs_gff

    - `args.direction` determines the frame scanning function:
        - "both" : scan_frame_both
        - "fw"   : scan_frame_fw

    - `args.fasta_output` determines which writer function to use when
      FASTA output is requested:
        - "nuc" : nucleotide output
        - "aa"  : amino acid output (requires codon table generation)

After configuring the functions, the script calls the selected ORF
prediction function with the provided input FASTA file and output base
name.
"""
if args.format == "both":
    scan_sequence = scan_sequence_both
elif args.format == "fasta":
    scan_sequence = scan_sequence_fasta
elif args.format == "gff":
    scan_sequence = scan_sequence_gff
    writer = writer_gff

if args.direction == "both":
    scan_frame = scan_frame_both
elif args.direction == "fw":
    scan_frame = scan_frame_fw

if args.fasta_output == "nuc" and args.format == "both":
    writer = writer_both_nuc
elif args.fasta_output == "nuc" and args.format == "fasta":
    writer = writer_fasta_nuc
elif args.fasta_output == "aa" and args.format == "both":
    codon_2_aa_fw,codon_2_aa_rv = const.CODON_2_AA_FW, const.CODON_2_AA_RV
    writer = writer_both_aa
elif args.fasta_output == "aa" and args.format == "fasta":
    codon_2_aa_fw,codon_2_aa_rv = const.CODON_2_AA_FW, const.CODON_2_AA_RV
    writer = writer_fasta_aa

"""
---------------------
ARGPARSE CUTOFF VALUE
---------------------
"""
"""
This section converts the user-defined cutoff value from the argument
parser into a nucleotide-based length threshold. 
The cutoff specified by the user represents the minimum
ORF length in amino acids.
"""

codon_length = 3
stop_codon_length = 3
cutoff_value = int(args.cutoff) * codon_length + stop_codon_length

"""
---------------------
ARGPARSE START CODONS
---------------------
"""
"""
This section processes the `--start_codons` argument from the command-line
parser and defines which start codons should be considered during ORF
detection. The canonical bacterial start codon ATG can optionally be
extended by including the alternative start codons GTG and/or TTG.

For reverse-strand scanning, the corresponding reverse-complement start
codons are also defined.
"""

if args.start_codons == "ATG":
    start_codons = ["ATG"]
    start_codons_rv = ["CAT"]
elif args.start_codons == "+GTG":
    print("Note that alternative start codons are often transcribed as methionin.")
    start_codons = ["ATG","GTG"]
    start_codons_rv = ["CAT","CAC"]
elif args.start_codons == "+TTG":
    print("Note that alternative start codons are often transcribed as methionin.")
    start_codons = ["ATG","TTG"]
    start_codons_rv = ["CAT","AAC"]
elif args.start_codons == "+GTG+TTG":
    print("Note that alternative start codons are often transcribed as methionin.")
    start_codons = ["ATG","GTG","TTG"]
    start_codons_rv = ["CAT","CAC","AAC"]

"""
------------
START SCRIPT
------------
"""

read_fasta(args.input, args.output)
