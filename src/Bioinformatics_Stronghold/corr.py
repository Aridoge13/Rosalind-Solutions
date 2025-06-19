# Error Correction in Reads
# URL: https://rosalind.info/problems/corr/

# Given: A collection of up to 1000 reads of equal length (at most 50 bp) in FASTA format. Some of these reads were generated with a single-nucleotide error. For each read s
# in the dataset, one of the following applies:
# 1. s was correctly sequenced and appears in the dataset at least twice(possibly as a reverse complement)
# 2. s is incorrect, it appears in the dataset exactly once, and its Hamming distance is 1 with respect to exactly one correct read in the dataset (or its reverse complement).

# Return: A list of all corrections in the form "[old read]->[new read]". (Each correction must be a single symbol substitution, and you may return the corrections in any order.)


from Bio import SeqIO
from collections import Counter
from Bio.Seq import Seq


def read_fasta(file_path):
    """Read FASTA file and return a list of reads."""
    return [str(record.seq) for record in SeqIO.parse(file_path, "fasta")]


def hamming_distance(s1, s2):
    """Calculate the Hamming distance between two strings."""
    if len(s1) != len(s2):
        raise ValueError("Strings must be of equal length")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))


def find_corrections(reads):

    # Count all reads and their reverse complements
    all_reads = reads + [str(Seq(r).reverse_complement()) for r in reads]
    read_counter = Counter(all_reads)

    # Correct reads: those that appear at least twice (including reverse complements)
    correct_reads = set(r for r in reads if read_counter[r] > 1 or read_counter[str(
        Seq(r).reverse_complement())] > 1)

    corrections = []
    for r in reads:
        if read_counter[r] == 1 and read_counter[str(Seq(r).reverse_complement())] == 1:
            # Try to correct to any correct read or its reverse complement
            for c in correct_reads:
                if hamming_distance(r, c) == 1:
                    corrections.append(f"{r}->{c}")
                    break
                rc = str(Seq(c).reverse_complement())
                if hamming_distance(r, rc) == 1:
                    corrections.append(f"{r}->{rc}")
                    break
    return corrections


if __name__ == "__main__":
    input_path = input("Enter the path to the FASTA file: ")
    reads = read_fasta(input_path)
    corrections = find_corrections(reads)
    output_path = input("Enter the path to save corrections: ")
    with open(output_path, 'w') as f:
        for correction in corrections:
            f.write(correction + '\n')
    print(f"Corrections saved to {output_path}")
