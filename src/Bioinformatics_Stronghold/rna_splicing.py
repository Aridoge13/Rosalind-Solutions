# RNA Splicing
# URL: https://rosalind.info/problems/splc/

# Given:  A DNA string s(of length at most 1 kbp) and a collection of substrings of s acting as introns. All strings are given in FASTA format.
# Return:  A protein string resulting from transcribing and translating the exons of s. (Note: Only one solution will exist for the dataset provided.)

# MVYIADKQHVASREAYGHMFKVCA
# MVYIADKQHVASREAYGHMFKVCA

from Bio import SeqIO
from Bio.Seq import Seq


def rna_splicing(fasta_file):
    # Read the DNA sequence and introns from the FASTA file
    with open(fasta_file, "r") as file:
        records = list(SeqIO.parse(file, "fasta"))
        dna_sequence = str(records[0].seq)
        introns = [str(record.seq) for record in records[1:]]

    # Remove introns from the DNA sequence
    for intron in introns:
        dna_sequence = dna_sequence.replace(intron, "")

    # Transcribe DNA to RNA
    rna_sequence = dna_sequence.replace("T", "U")

    # Translate RNA to protein using Biopython
    protein = str(Seq(rna_sequence).translate(to_stop=True))

    return protein


if __name__ == "__main__":
    fasta_file = input("Enter the path to the FASTA file: ").strip()
    protein_sequence = rna_splicing(fasta_file)
    output_file = input("Enter the path to the output file: ").strip()
    with open(output_file, "w") as file:
        file.write(protein_sequence)
    print(f"Protein sequence written to {output_file}")
