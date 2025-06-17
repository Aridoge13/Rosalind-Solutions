# Open Reading Frame
# URL: https://rosalind.info/problems/orf/

# Given:  A DNA string s of length at most 1 kbp in FASTA format
# Result: Every distinct candidate protein string that can be translated from ORFs of s. Strings can be returned in any order.

from Bio.Seq import Seq
from Bio import SeqIO

STOP_CODONS = {"TAA", "TAG", "TGA"}


def find_orfs_in_seq(seq: Seq) -> set:
    proteins = set()
    seq_len = len(seq)
    for frame in range(3):
        i = frame
        while i < seq_len - 2:
            codon = str(seq[i:i+3])
            if codon == "ATG":
                for j in range(i + 3, seq_len - 2, 3):
                    stop_codon = str(seq[j:j+3])
                    if stop_codon in STOP_CODONS:
                        orf_dna = seq[i:j+3]
                        protein = orf_dna.translate(to_stop=True)
                        proteins.add(str(protein))
                        break  # stop after the first in-frame stop codon
            i += 3  # Move to next codon
    return proteins


def find_all_orfs(dna: str) -> set:
    dna_seq = Seq(dna)
    all_proteins = set()
    strands = [dna_seq, dna_seq.reverse_complement()]
    for strand in strands:
        all_proteins.update(find_orfs_in_seq(strand))
    return all_proteins


def read_fasta(filepath: str) -> str:
    record = SeqIO.read(filepath, "fasta")
    return str(record.seq)


if __name__ == "__main__":
    fasta_file = input("Enter the path to the FASTA file: ").strip()
    dna_seq = read_fasta(fasta_file)
    proteins = find_all_orfs(dna_seq)
    for p in proteins:
        print(p)
