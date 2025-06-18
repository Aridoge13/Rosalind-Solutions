# Locating Restriction Sites
# URL: https://rosalind.info/problems/revp/

# Given: A DNA string of length at most 1 kbp in FASTA format.
# Return: The position and length of every reverse palindrome in the string having length between 4 and 12. You may return these pairs in any order.

from Bio import SeqIO
from Bio.Seq import Seq


def find_reverse_palindromes(dna, min_length=4, max_length=12):
    palindromes = []
    n = len(dna)

    for length in range(min_length, max_length + 1):
        for i in range(n - length + 1):
            substring = dna[i:i + length]
            # Check if the substring is a palindrome
            if substring == str(Seq(substring).reverse_complement()):
                # Store position (1-based index) and length
                palindromes.append((i + 1, length))

    return palindromes


if __name__ == "__main__":
    # Read the DNA sequence from a FASTA file
    fasta_file = input("Enter the path to the FASTA file: ").strip()
    with open(fasta_file, "r") as file:
        record = SeqIO.read(file, "fasta")
        dna_sequence = str(record.seq)

    # Find reverse palindromes
    results = find_reverse_palindromes(dna_sequence)

    # Print results
    for position, length in results:
        print(position, length)
