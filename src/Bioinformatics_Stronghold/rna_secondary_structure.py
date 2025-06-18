# Perfect Matchings and RNA Secondary Structures
# URL: https://rosalind.info/problems/pmch/

# Given: An RNA string s of length at most 80 bp having the same number of occurrences of 'A' as 'U' and the same number of occurrences of 'C' as 'G'.
# Return: The total possible number of perfect matchings of basepair edges in the bonding graph of s.

import math


def count_perfect_matchings(s):
    a = s.count('A')
    u = s.count('U')
    c = s.count('C')
    g = s.count('G')
    # The problem guarantees that a == u and c == g
    return math.factorial(a) * math.factorial(c)


if __name__ == "__main__":
    rna_string = input("Enter path to RNA string file: ").strip()
    with open(rna_string, 'r') as file:
        rna_string = ''.join(line.strip()
                             for line in file if not line.startswith('>'))
    output_path = input("Enter path to output file: ").strip()
    result = count_perfect_matchings(rna_string)
    with open(output_path, 'w') as file:
        file.write(str(result) + '\n')
