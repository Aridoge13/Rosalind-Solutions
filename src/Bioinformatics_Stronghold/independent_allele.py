# Independent Alleles
# URL: https://rosalind.info/problems/lia/

# Given: Two positive integers k(k≤7) and N(N≤2k). In this problem, we begin with Tom, who in the 0th generation has genotype Aa Bb. Tom has two children in the 1st generation, each of whom has two children, and so on. Each organism always mates with an organism having genotype Aa Bb.
# Return: The probability that at least N Aa Bb organisms will belong to the k-th generation of Tom's family tree (don't count the Aa Bb mates at each level). Assume that Mendel's second law holds for the factors.


def independent_allele(k, n):
    total_organisms = 2 ** k  # Total number of organisms in the k-th generation

    if n > total_organisms:
        return 0.0

    from math import comb

    probability = sum(
        comb(total_organisms, i) * (0.25 ** i) *
        (0.75 ** (total_organisms - i))
        for i in range(n, total_organisms + 1)
    )

    return probability  # Probability of at least N AaBb organisms in the k-th generation


if __name__ == "__main__":
    input_data = input("Enter k and N separated by space: ")
    k, n = map(int, input_data.split())
    result = independent_allele(k, n)
    print(
        f"The probability of at least {n} AaBb organisms in the {k}-th generation is: {result:.5f}")
