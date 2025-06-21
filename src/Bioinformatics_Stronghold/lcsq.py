# Finding a Shared Spliced Motif
# URL: https://rosalind.info/problems/lcsq/

# Given: Two DNA strings s and t (each of length at most 1 kbp) in FASTA format.
# Return: A longest common subsequence of s and t. (If more than one solution exists, you may return any one.)

from Bio import SeqIO
from typing import List


def read_fasta(file_path: str) -> List[str]:
    """Reads a FASTA file and returns the sequences as a list of strings."""
    sequences = []
    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, "fasta"):
            sequences.append(str(record.seq))
    return sequences


# Sample Input:
# >Rosalind_23
# AACCTTGG
# >Rosalind_64
# ACACTGTGA

def longest_common_subsequence(s: str, t: str) -> str:
    n, m = len(s), len(t)
    dp = [[""] * (m + 1) for _ in range(n + 1)]
    for i in range(n):
        for j in range(m):
            if s[i] == t[j]:
                dp[i + 1][j + 1] = dp[i][j] + s[i]
            else:
                dp[i + 1][j + 1] = max(dp[i][j + 1], dp[i + 1][j], key=len)
    return dp[n][m]


if __name__ == "__main__":
    input_file = input("Enter the path to the FASTA file: ").strip()
    sequences = read_fasta(input_file)

    if len(sequences) != 2:
        raise ValueError("Expected exactly two sequences in the FASTA file.")

    s, t = sequences
    lcs_result = longest_common_subsequence(s, t)

    output_path = input("Enter the path to save the output: ").strip()
    with open(output_path, 'w') as output_file:
        output_file.write(lcs_result + '\n')
