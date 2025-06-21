# k-Mer Composition
# URL: http://rosalind.info/problems/kmer/

# Given: A DNA string s in FASTA format (having length at most 1000 bp).
# Return: The 4-mer composition of s.

from itertools import product


def kmer_composition(s):
    # Generate all possible 4-mers in lex order
    bases = ['A', 'C', 'G', 'T']
    all_kmers = [''.join(p) for p in product(bases, repeat=4)]
    kmer_count = {kmer: 0 for kmer in all_kmers}

    for i in range(len(s) - 3):
        kmer = s[i:i + 4]
        if kmer in kmer_count:
            kmer_count[kmer] += 1

    return [kmer_count[kmer] for kmer in all_kmers]


def read_fasta(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return ''.join(line.strip() for line in lines if not line.startswith('>'))


if __name__ == "__main__":
    input_file = input("Enter the path to the FASTA file: ").strip()
    dna_string = read_fasta(input_file)
    kmer_counts = kmer_composition(dna_string)
    output_path = input("Enter the path to save the output: ").strip()
    with open(output_path, 'w') as output_file:
        output_file.write(' '.join(map(str, kmer_counts)) + '\n')
