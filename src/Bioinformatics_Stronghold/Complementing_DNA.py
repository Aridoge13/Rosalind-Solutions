# Complementing a Strand of DNA

# Problem: In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C' and 'G'.
# The reverse complement of a DNA string s is the string sc formed by reversing the symbols of s, then taking the complement of each symbol(e.g., the reverse complement of "GTCA" is "TGAC").

# Given: A DNA string s of length at most 1000 bp.
# Return: The reverse complement sc of s.


def reverse_complement(dna_string):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reversed_string = dna_string[::-1]
    reverse_complement_string = ''.join(
        complement[base] for base in reversed_string)
    return reverse_complement_string


if __name__ == "__main__":
    input_path = input("Enter the path to your input file here: ")
    with open(input_path, 'r') as file:
        input_dna = file.read().strip()

    output_path = input("Enter the path to your output file here: ")
    with open(output_path, 'w') as file:
        file.write(reverse_complement(input_dna))
    result = reverse_complement(input_dna)
    print("Reverse complement:", result)
