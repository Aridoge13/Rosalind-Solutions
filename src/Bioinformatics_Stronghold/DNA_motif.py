# Finding a motif in DNA
# URL: https://rosalind.info/problems/subs/
# Combing through the haystack
# Given: Two DNA strings s and t (each of length at most 1 kbp).
# Return: All locations of t as a substring of s (if any such locations exist).

from Bio import SeqIO


def find_motif_locations(s, t):
    locations = []
    start = 0

    while True:
        start = s.find(t, start)
        if start == -1:
            break
        locations.append(start + 1)  # Convert to 1-based index
        start += 1  # Move to the next character after the current match

    return locations


if __name__ == "__main__":
    # Read the input strings from a txt file
    input_path = input("Enter the path to the input file: ")
    with open(input_path, 'r') as file:
        lines = file.readlines()
        s = lines[0].strip()  # First line is the DNA string s
        t = lines[1].strip()  # Second line is the motif t
    locations = find_motif_locations(s, t)
    if locations:
        print(" ".join(map(str, locations)))
    else:
        print("No motif found.")
