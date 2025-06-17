# Consensus and Profile
# URL: https://rosalind.info/problems/cons/
# Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.
# Return: A consensus string and profile matrix for the collection. (If several possible consensus strings exist, then you may return any one of them.)

from Bio import SeqIO


def consensus_profile(fasta_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    length = len(sequences[0].seq)

    # Initialize profile matrix
    profile = {nucleotide: [0] * length for nucleotide in "ACGT"}

    # Build the profile matrix
    for seq in sequences:
        for i, nucleotide in enumerate(seq.seq):
            profile[nucleotide][i] += 1

    # Build the consensus string
    consensus = ''.join(
        max("ACGT", key=lambda x: profile[x][i]) for i in range(length))

    return consensus, profile


if __name__ == "__main__":
    fasta_file = input("Enter the path to the FASTA file: ")
    consensus, profile = consensus_profile(fasta_file)
    print("Consensus String:", consensus)   # Print the consensus string
    output_path = input("Enter the path to save the profile matrix: ")
    with open(output_path, "w") as f:
        for nucleotide in "ACGT":
            f.write(
                f"{nucleotide}: {' '.join(map(str, profile[nucleotide]))}\n")
    print(f"Profile matrix saved to {output_path}")


# The consensus string will be printed to the terminal, and the profile matrix will be saved to a specified file.
# Example usage:
# python consensus_profile.py
# Enter the path to the FASTA file: input.fasta
# Consensus String: ACGT...
# Enter the path to save the profile matrix: output.txt
# Profile matrix saved to output.txt
