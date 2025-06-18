# Enumerating Gene Orders
# URL: https://rosalind.info/problems/perm/

# Given: A positive integer n (n ≤ 7).
# Return: The total number of permutations of length n, followed by a list of all permutations (in any order).


from itertools import permutations


def enumerate_gene_orders(n):
    # Generate all permutations of the numbers 1 to n
    perm = permutations(range(1, n + 1))
    # Convert permutations from tuples to lists
    perm_list = [list(p) for p in perm]
    # Return the total number of permutations and the list of permutations
    return len(perm_list), perm_list


if __name__ == "__main__":
    n = int(input("Enter your input here:").strip())
    total_permutations, permutations_list = enumerate_gene_orders(n)

    # Print the total number of permutations
    print(total_permutations)
    # Print each permutation on a new line
    for perm in permutations_list:
        print(' '.join(map(str, perm)))
