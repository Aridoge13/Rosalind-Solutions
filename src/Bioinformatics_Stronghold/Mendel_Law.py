# Mendel's First Law
# URL: http://rosalind.info/problems/iprb/
# Mendel's First Law states that the probability of an offspring being dominant or recessive can be calculated based on the genotypes of the parents.

def mendel_law(k, m, n):
    """
    Calculate the probability of producing a dominant phenotype offspring
    given the number of homozygous dominant (k), heterozygous (m), and homozygous recessive (n) individuals.
    :param k: Number of homozygous dominant individuals (AA)
    :param m: Number of heterozygous individuals (Aa)
    :param n: Number of homozygous recessive individuals (aa)
    :return: Probability of producing a dominant phenotype offspring
    """
    total = k + m + n
    if total == 0:
        return 0.0

    total_pairs = total * (total - 1)

    # Each pairing and its probability of producing a dominant offspring
    prob = 0
    prob += k * (k - 1) * 1.0          # AA x AA
    prob += k * m * 2 * 1.0            # AA x Aa and Aa x AA
    prob += k * n * 2 * 1.0            # AA x aa and aa x AA
    prob += m * (m - 1) * 0.75         # Aa x Aa
    prob += m * n * 2 * 0.5            # Aa x aa and aa x Aa
    # aa x aa gives 0, so not included

    return prob / total_pairs


# Provide the path to the input file containing k, m, n values
if __name__ == "__main__":
    input_path = input("Enter the path to the input file: ")
    with open(input_path, 'r') as file:
        k, m, n = map(int, file.readline().strip().split())
    result = mendel_law(k, m, n)
    print(result)


# Example Input:
# 2 2 2
# Example Output:
# 0.7833333333333333
