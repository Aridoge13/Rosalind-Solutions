# Partial Permutations
# URL: https://rosalind.info/problems/pper/

# Given: Positive integers n and k such that 100 ≥ n > 0 and 10 ≥ k > 0.
# Return: The total number of partial permutations P(n, k), modulo 1,000,000.

def partial_permutations(n, k):
    result = 1
    for i in range(n, n - k, -1):
        result = (result * i) % 1000000
    return result


if __name__ == "__main__":
    n, k = map(int, input().strip().split())
    print(partial_permutations(n, k))
