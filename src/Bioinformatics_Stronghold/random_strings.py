# Introduction to Random Strings
# URL: https://rosalind.info/problems/prob/

# Given: A DNA string s of length at most 1000 bp and an array A containing at most 20 numbers between 0 and 1.
# Return: An array B having the same length as A in which B[k] represents the common logarithm of the probability that a random string constructed with the GC-content found in A[k] will match s exactly.


from Bio import SeqIO
from math import log10
from typing import List

# Example input:
# s = ACGATACAA
# A = 0.129 0.287 0.423 0.476 0.641 0.742 0.783
# Output for the example input:
# -5.737 -5.217 -5.263 -5.360 -5.958 -6.628 -7.009


def calculate_log_probabilities(s: str, A: List[float]) -> List[float]:
    """
    Calculate the common logarithm of the probability that a random string constructed with the GC-content found in A will match s exactly.

    :param s: DNA string
    :param A: List of GC-content values
    :return: List of log probabilities
    """
    length = len(s)
    probabilities = []

    for gc_content in A:
        gc_count = sum(1 for base in s if base in 'GC')
        at_count = length - gc_count

        p_gc = gc_content / 2
        p_at = (1 - gc_content) / 2

        probability = (p_gc ** gc_count) * (p_at ** at_count)
        log_probability = log10(probability)

        probabilities.append(log_probability)

    return probabilities


def read_fasta(file_path: str) -> str:
    """
    Read a FASTA file and return the DNA string.

    :param file_path: Path to the FASTA file
    :return: DNA string
    """
    with open(file_path, 'r') as file:
        record = SeqIO.read(file, 'fasta')
        return str(record.seq)


if __name__ == "__main__":
    import sys
    lines = sys.stdin.read().strip().splitlines()
    # Read DNA string from FASTA
    s = ''.join(line.strip() for line in lines if not line.startswith(
        '>') and set(line.strip()) <= set('ACGT'))
    # Read GC-content array (last line)
    A = list(map(float, lines[-1].strip().split()))
    log_probabilities = calculate_log_probabilities(s, A)
    print(" ".join(f"{prob:.3f}" for prob in log_probabilities))
