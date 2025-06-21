# Speeding Up Motif Finding
# URL: https://rosalind.info/problems/kmp/

# Given: A DNA string s (of length at most 100 kbp) in FASTA format.
# Return: The failure array for s.

from Bio import SeqIO
from typing import List


def read_fasta(file_path: str) -> str:
    """Reads a FASTA file and returns the sequence as a string."""
    with open(file_path, 'r') as file:
        record = SeqIO.read(file, "fasta")
        return str(record.seq)


def failure_array(s: str) -> List[int]:
    """Computes the failure array for the given string s."""
    n = len(s)
    f = [0] * n
    j = 0  # length of previous longest prefix suffix

    for i in range(1, n):
        while j > 0 and s[i] != s[j]:
            j = f[j - 1]
        if s[i] == s[j]:
            j += 1
        f[i] = j

    return f


if __name__ == "__main__":
    input_path = input("Enter the path to the FASTA file: ")
    s = read_fasta(input_path)
    result = failure_array(s)
    output_path = input("Enter the path to save the output: ")
    with open(output_path, 'w') as output_file:
        output_file.write(' '.join(map(str, result)))
