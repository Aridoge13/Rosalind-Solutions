# Introduction to the Bioinformatics Armory
# URL: https://rosalind.info/problems/ini/

# Given: A DNA string s of length at most 1000 bp.
# Return: Four integers (separated by spaces) representing the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s.


def count_nucleotide_bases(dna_string):
    bases = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    for nucleotide in dna_string:
        if nucleotide in bases:
            bases[nucleotide] += 1
    return [bases['A'], bases['C'], bases['G'], bases['T']]


if __name__ == "__main__":
    dna_string = input("Please input the DNA sequence: ").strip()
    result = count_nucleotide_bases(dna_string)
    output_path = input("Please provide the path for the output file: ")
    with open(output_path, 'w') as file:
        file.write(' '.join(map(str, result)) + '\n')
