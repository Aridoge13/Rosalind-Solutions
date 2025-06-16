# Problem 1: Counting DNA Nucleotides
# Given: A DNA string s of length at most 1000 nt.
# Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s

import sys


def count_nucleotides(dna_string):
    # Initialize counts for each nucleotide
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

    # Iterate through each character in the DNA string
    for nucleotide in dna_string:
        if nucleotide in counts:
            counts[nucleotide] += 1

    # Return the counts as a list
    return [counts['A'], counts['C'], counts['G'], counts['T']]


# Input the path to the DNA file

if __name__ == "__main__":
    input_file = input("Enter the path to the DNA file: ").strip()

    try:
        # Read the DNA string from the file
        with open(input_file, 'r') as file:
            dna_string = file.read().strip()

        # Count nucleotides
        result = count_nucleotides(dna_string)

        # Print the result
        print(" ".join(map(str, result)))

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
