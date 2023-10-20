import argparse
import sys
import os
import csv
import re
from re import Pattern
from pathlib import Path
from typing import List, Union, Optional
import textwrap


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='genome_file', type=isfile, required=True, 
                        help="Complete genome file in fasta format")
    parser.add_argument('-g', dest='min_gene_len', type=int, 
                        default=50, help="Minimum gene length to consider (default 50).")
    parser.add_argument('-s', dest='max_shine_dalgarno_distance', type=int, 
                        default=16, help="Maximum distance from start codon "
                        "where to look for a Shine-Dalgarno motif (default 16).")
    parser.add_argument('-d', dest='min_gap', type=int, default=40,
                        help="Minimum gap between two genes - shine box not included (default 40).")
    parser.add_argument('-p', dest='predicted_genes_file', type=Path, 
                        default=Path("predict_genes.csv"),
                        help="Tabular file giving position of predicted genes")
    parser.add_argument('-o', dest='fasta_file', type=Path,
                        default=Path("genes.fna"),
                        help="Fasta file giving sequence of predicted genes")
    return parser.parse_args()


def read_fasta(fasta_file: Path) -> str:
    """Extract genome sequence from fasta files.

    :param fasta_file: (Path) Path to the fasta file.
    :return: (str) Sequence from the genome. 
    """
    with open(fasta_file, "r") as fasta_input:
        #Skip header
        fasta_input.readline()
        #Get genome sequence
        sequence = ""
        for ligne in fasta_input:
            sequence += ligne.strip().upper()
    return sequence


def find_start(start_regex: Pattern, sequence: str, start: int, stop: int) -> Union[int, None]:
    """Find next start codon before a end position.

    :param start_regexp: A regex object that identifies a start codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :param stop: (int) Stop position of the research
    :return: (int) If exist, position of the start codon. Otherwise None. 
    """
    match_start = start_regex.search(sequence, pos = start, endpos = stop)
    if match_start:
        #Return start position if macth found
        return match_start.start(0)
    return None

def find_stop(stop_regex: Pattern, sequence: str, start: int) -> Union[int, None]:
    """Find next stop codon that should be in the same reading phase as the start.

    :param stop_regexp: A regex object that identifies a stop codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :return: (int) If exist, position of the stop codon. Otherwise None. 
    """
    #Start position is provided by find start
    #Redefine start in order to skip start codon
    #Looking stop codon patterns from start
    match_stops = stop_regex.finditer(sequence, pos = start)
    for match_stop in match_stops:
        start_pos = match_stop.start(0) 
        if (start_pos - start) % 3 == 0:
            #Return first macth found
            #Codon has a size of 3. We need to use modulo 3 to detect a codon
            return start_pos
    return None


def has_shine_dalgarno(shine_regex: Pattern, sequence: str, start: int, max_shine_dalgarno_distance: int) -> bool:
    """Find a shine dalgarno motif before the start codon

    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Position of the start in the genome
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine dalgarno to the start position
    :return: (boolean) true -> has a shine dalgarno upstream to the gene, false -> no
    """
    #start position is provided by find start
    #Define start position and end position to look
    end_pos = start - 6 #start codon position
    start_pos = start - max_shine_dalgarno_distance #max distance before start
    if start_pos < 0:
        return False
    #Not need end_pos + 1 beacause it correspond to codon start first position
    match_shine = shine_regex.search(sequence, pos = start_pos, endpos = end_pos)
    if match_shine:
        return True
    return False
 
def predict_genes(sequence: str, start_regex: Pattern, stop_regex: Pattern, shine_regex: Pattern, 
                  min_gene_len: int, max_shine_dalgarno_distance: int, min_gap: int) -> List[List[int]]:
    """Predict most probable genes

    :param sequence: (str) Sequence from the genome.
    :param start_regexp: A regex object that identifies a start codon.
    :param stop_regexp: A regex object that identifies a stop codon.
    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param min_gene_len: (int) Minimum gene length.
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine dalgarno to the start position.
    :param min_gap: (int) Minimum distance between two genes.
    :return: (list) List of [start, stop] position of each predicted genes.
    """
    #Create empty gene_pos list
    gene_pos = []
    #Define current opsition variable
    current_pos = 0
    #Brwose the genome
    while len(sequence) - current_pos >= min_gap:
        #Looking for start codon position
        current_pos = find_start(start_regex, sequence, current_pos, len(sequence))
        if current_pos:
            #Looking for stop codon position
            stop = find_stop(stop_regex, sequence, current_pos)
            if stop:
                gene_length = stop - current_pos #define gene length
                if gene_length >= min_gene_len:
                    #Looking shine dalgaro pattern before start codon
                    match_shine = has_shine_dalgarno(shine_regex, sequence, current_pos, max_shine_dalgarno_distance)
                    if match_shine:
                        #Must identify a gene
                        #Append gene position to the list
                        gene_pos.append([current_pos + 1, stop + 3]) #+1 and +3 for teh true position
                        #Update start position
                        current_pos = stop + 3 + min_gap
                    else:
                        #May be we don't identify a gene
                        current_pos += 1
                else:
                    current_pos += 1
            else:
                current_pos += 1
    return gene_pos

def write_genes_pos(predicted_genes_file: Path, probable_genes: List[List[int]]) -> None:
    """Write list of gene positions.

    :param predicted_genes_file: (Path) Output file of gene positions.
    :param probable_genes: List of [start, stop] position of each predicted genes.
    """
    try:
        with predicted_genes_file.open("wt") as predict_genes:
            predict_genes_writer = csv.writer(predict_genes, delimiter=",")
            predict_genes_writer.writerow(["Start", "Stop"])
            predict_genes_writer.writerows(probable_genes)
    except IOError:
        sys.exit("Error cannot open {}".format(predicted_genes_file))


def write_genes(fasta_file: Path, sequence: str, probable_genes: List[List[int]], sequence_rc: str, 
                probable_genes_comp: List[List[int]]):
    """Write gene sequence in fasta format

    :param fasta_file: (Path) Output fasta file.
    :param sequence: (str) Sequence of genome file in 5'->3'.
    :param probable_genes: (list) List of [start, stop] position of each predicted genes in 5'->3'.
    :param sequence_rc: (str) Sequence of genome file in 3' -> 5'.
    :param probable_genes_comp: (list)List of [start, stop] position of each predicted genes in 3' -> 5'.
    """
    try:
        with open(fasta_file, "wt") as fasta:
            for i,gene_pos in enumerate(probable_genes):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                    i+1, os.linesep, 
                    textwrap.fill(sequence[gene_pos[0]-1:gene_pos[1]])))
            i = i+1
            for j,gene_pos in enumerate(probable_genes_comp):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                            i+1+j, os.linesep,
                            textwrap.fill(sequence_rc[gene_pos[0]-1:gene_pos[1]])))
    except IOError:
        sys.exit("Error cannot open {}".format(fasta_file))


def reverse_complement(sequence: str) -> str:
    """Get the reverse complement

    :param sequence: (str) DNA Sequence.
    :return: (str) Reverse complemented sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in sequence[::-1]])


#==============================================================
# Main program
#==============================================================
def main() -> None: # pragma: no cover
    """
    Main program function
    """
    # Gene detection over genome involves to consider a thymine instead of
    # an uracile that we would find on the expressed RNA
    #start_codons = ['TTG', 'CTG', 'ATT', 'ATG', 'GTG']
    #stop_codons = ['TAA', 'TAG', 'TGA']
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')
    # Shine AGGAGGUAA
    #AGGA ou GGAGG 
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')
    # Arguments
    args = get_arguments()
    # Get genome sequence
    print("Get genome sequence from FASTA file")
    sequence = read_fasta(args.genome_file)
    # Predict gene position in sequence
    print("Predict gene position in 5'-> 3' direction")
    gene_list = predict_genes(sequence, start_regex, stop_regex, shine_regex, 
    args.min_gene_len, args.max_shine_dalgarno_distance, args.min_gap)
    # Call reverse complement
    print("Get reverse complement sequence")
    sequence_rc = reverse_complement(sequence)

    print("Predict gene position in 3'-> 5' direction")
    gene_list_reverse = predict_genes(sequence, start_regex, stop_regex, shine_regex,
    args.min_gene_len, args.max_shine_dalgarno_distance, args.min_gap)

    #Modify gene_list_reverse
    print("Combine gene position results")
    for i in range(len(gene_list_reverse)):
        for j in range(len(gene_list_reverse[i])):
            gene_list_reverse[i][j] = len(sequence) - gene_list_reverse[i][j] 
    #Combine list
    gene_list_total = gene_list + gene_list_reverse
    
    # Write results
    print("Write results")
    write_genes_pos(args.predicted_genes_file, gene_list_total)
    write_genes(args.fasta_file, sequence, gene_list, sequence_rc, gene_list_reverse)


if __name__ == '__main__':
    main()
