# Inferring mRNA from a Protein Sequence
# URL: https://rosalind.info/problems/mrna/

# Given: A protein string of length at most 1000 aa.
# Result: The total number of different RNA strings from which the protein could have been translated, modulo 1,000,000. (Don't neglect the importance of the stop codon in protein translation.)

from Bio.Data import CodonTable
from collections import Counter


def protein_mrna(protein: str) -> int:
    """
    Calculate the number of different mRNA strings that can encode a given protein sequence.
    """
    table = CodonTable.unambiguous_rna_by_name["Standard"]
    codon_count = Counter(table.forward_table.values())
    codon_count['*'] = len(table.stop_codons)
    mod = 1_000_000
    total = 1

    for aa in protein:
        total = (total * codon_count[aa]) % mod

    # Multiply by stop codon count at the end
    total = (total * codon_count['*']) % mod
    return total


if __name__ == "__main__":
    protein_file = input(
        "Enter the path to the file containing the protein sequence: ").strip()
    with open(protein_file, 'r') as file:
        protein_sequence = file.read().strip()
    result = protein_mrna(protein_sequence)
    print(result)
