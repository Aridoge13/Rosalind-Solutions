# Counting Phylogenetic Ancestors
# URL: https://rosalind.info/problems/inod/

# Given: A positive integer n (3 ≤ n ≤ 1000).
# Return: The number of internal nodes of any unrooted binary tree having n leaves

def count_internal_nodes(n):
    # The number of internal nodes in an unrooted binary tree with n leaves is n - 2
    return n - 2


if __name__ == "__main__":
    n = int(input().strip())
    result = count_internal_nodes(n)
    print(result)
