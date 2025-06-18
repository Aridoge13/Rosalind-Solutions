# Genome Assembly as Shortest Superstring
# URL: https://rosalind.info/problems/long/

# Given: At most 50 DNA strings of approximately equal length, not exceeding 1 kbp, in FASTA format (which represent reads deriving from the same strand of a single linear chromosome).
# The dataset is guaranteed to satisfy the following condition: there exists a unique way to reconstruct the entire chromosome from these reads by gluing together pairs of reads that overlap by more than half their length.
# Return: A  shortest superstring containing all the given strings (thus corresponding to a reconstructed chromosome).

from Bio import SeqIO


def read_fasta(file_path):
    """Read a FASTA file and return a list of sequences."""
    sequences = []
    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            sequences.append(str(record.seq))
    return sequences


def overlap(s1, s2):
    """Return the length of the overlap between two strings."""
    max_overlap = 0
    min_length = min(len(s1), len(s2))
    for i in range(1, min_length + 1):
        if s1[-i:] == s2[:i]:
            max_overlap = i
    return max_overlap


def merge_strings(s1, s2):
    """Merge two strings with the maximum overlap."""
    ov = overlap(s1, s2)
    return s1 + s2[ov:]


def shortest_superstring(sequences):
    while len(sequences) > 1:
        max_overlap = -1
        best_pair = (0, 0)
        best_direction = 0  # 0: i->j, 1: j->i
        for i in range(len(sequences)):
            for j in range(len(sequences)):
                if i == j:
                    continue
                ov = overlap(sequences[i], sequences[j])
                if ov > max_overlap:
                    max_overlap = ov
                    best_pair = (i, j)
                    best_direction = 0
                ov_rev = overlap(sequences[j], sequences[i])
                if ov_rev > max_overlap:
                    max_overlap = ov_rev
                    best_pair = (j, i)
                    best_direction = 1

        if max_overlap == -1:
            break

        i, j = best_pair
        merged_string = merge_strings(sequences[i], sequences[j])
        # Remove both sequences and add the merged one
        sequences.pop(max(i, j))
        sequences.pop(min(i, j))
        sequences.append(merged_string)

    return sequences[0]


if __name__ == "__main__":
    input_file = input("Enter the path to the FASTA file: ")
    sequences = read_fasta(input_file)
    superstring = shortest_superstring(sequences)
    output_file = input("Enter the path to save the output: ")
    with open(output_file, 'w') as file:
        file.write(superstring + '\n')
