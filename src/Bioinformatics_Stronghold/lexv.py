# Ordering Strings of Varying Length Lexicographically
# URL: https://rosalind.info/problems/lexv/

# Given: A permutation of at most 12 symbols defining an ordered alphabet Z and a positive integer n (n≤4).
# Return: All strings of length n formed from Z ordered lexicographically. (Note: As in “Enumerating k-mers Lexicographically”, alphabet order is based on the order in which the symbols are given.)

from itertools import product


def lexv(alphabet, n):
    result = []
    for length in range(1, n + 1):
        for p in product(alphabet, repeat=length):
            result.append(''.join(p))

    # Sort using custom alphabet order
    order = {char: i for i, char in enumerate(alphabet)}
    result.sort(key=lambda word: [order[c] for c in word])
    return result


if __name__ == "__main__":
    # Example input:
    # A B C D
    # 2
    alphabet = input("Enter alphabet (space-separated): ").strip().split()
    n = int(input("Enter maximum length: ").strip())

    strings = lexv(alphabet, n)
    output_path = input("Enter the path to the output file: ")
    with open(output_path, 'w') as output_file:
        output_file.write('\n'.join(strings) + '\n')
