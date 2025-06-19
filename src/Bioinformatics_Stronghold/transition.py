# Transitions and Transversions
# URL: https://rosalind.info/problems/tran/

# Given: Two DNA strings s1 and s2 of equal length (at most 1 kbp).
# Return: The transition/transversion ratio R(s1,s2).

from Bio import SeqIO


def transition_transversion_ratio(s1, s2):
    transitions = 0
    transversions = 0

    for a, b in zip(s1, s2):
        if a != b:
            if (a in 'AG' and b in 'AG') or (a in 'CT' and b in 'CT'):
                transitions += 1
            else:
                transversions += 1

    if transitions + transversions == 0:
        return 0

    return transitions / transversions


def read_fasta(file_path):
    """
    Reads a FASTA file and returns a list of sequences.
    """
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    return sequences


if __name__ == "__main__":
    input_path = input("Enter the path to the input file: ")
    sequences = read_fasta(input_path)
    s1 = sequences[0]
    s2 = sequences[1]
    ratio = transition_transversion_ratio(s1, s2)
    print(f"{ratio:f}")
