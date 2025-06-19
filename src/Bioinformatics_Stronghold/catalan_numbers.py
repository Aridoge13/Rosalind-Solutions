# Catalan Numbers and RNA Secondary Structures
# URL: https://rosalind.info/problems/cat/

# Given: An RNA string s having the same number of occurrences of 'A' as 'U' and the same number of occurrences of 'C' as 'G'. The length of the string is at most 300 bp.
# Return: The total number of noncrossing perfect matchings of basepair edges in the bonding graph of s, modulo 1,000,000.

# HINT: Write a function that counts Catalan numbers via dynamic programming. How can we modify this function to apply to our given problem?

# sample input:
# > >Rosalind_57
# AUAU
# Output: 2

def count_matchings(rna, memo={}):
    if rna in memo:
        return memo[rna]
    if len(rna) == 0:
        return 1
    total = 0
    for i in range(1, len(rna), 2):
        if (rna[0], rna[i]) in [('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')]:
            left = count_matchings(rna[1:i], memo)
            right = count_matchings(rna[i+1:], memo)
            total += left * right
    memo[rna] = total % 1000000
    return memo[rna]


if __name__ == "__main__":
    input_path = input("Enter the path to the input file: ")
    with open(input_path, 'r') as file:
        lines = [line.strip() for line in file if not line.startswith('>')]
        rna_string = ''.join(lines)
    print(count_matchings(rna_string))
