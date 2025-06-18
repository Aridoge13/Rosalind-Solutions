# Enumerating K-mers Lexicographically
# URL: https://rosalind.info/problems/lexf/

# Given: A collection of at most 10 symbols defining an ordered alphabet, and a positive integer n (n ≤ 10).
# Return: All strings of length n that can be formed from the alphabet, ordered lexicographically (use the standard order of symbols in the English alphabet).

from itertools import product


def enumerate_kmers(alphabet, n):
    """
    Generate all strings of length n from the given alphabet, ordered lexicographically.

    :param alphabet: A list of symbols defining the ordered alphabet.
    :param n: The length of the strings to generate.
    :return: A list of strings of length n, ordered lexicographically.
    """
    return [''.join(p) for p in product(alphabet, repeat=n)]


def main():
    # Example input
    alphabet = ['A', 'C', 'G', 'T']
    n = 2

    # Generate k-mers
    kmers = enumerate_kmers(alphabet, n)

    # Print the result
    for kmer in kmers:
        print(kmer)


if __name__ == "__main__":
    input_alphabet = input("Enter the alphabet (space-separated): ").split()
    input_n = int(input("Enter the length of the strings (n): "))
    output_path = input("Enter the output file path: ").strip()
    with open(output_path, 'w') as file:
        for kmer in enumerate_kmers(input_alphabet, input_n):
            file.write(kmer + '\n')
