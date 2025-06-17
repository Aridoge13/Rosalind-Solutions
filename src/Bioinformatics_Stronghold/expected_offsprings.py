# Calculate expected offsprings
# URL: https://rosalind.info/problems/iev/
# Given: Six nonnegative integers, each of which does not exceed 20,000. The integers correspond to the number of couples in a population possessing each genotype pairing for a given factor. In order, the six given integers represent the number of couples having the following genotypes:
# AA-AA, AA-Aa, AA-aa, Aa-Aa, Aa-aa, aa-aa.
# Return: The expected number of offspring displaying the dominant phenotype in the next generation, under the assumption that every couple has exactly two offspring.

def expected_offsprings(genotypes):
    # Unpack the genotypes
    AA_AA, AA_Aa, AA_aa, Aa_Aa, Aa_aa, aa_aa = genotypes

    # Calculate the expected number of dominant phenotype offspring
    expected_dominant = (
        2 * AA_AA * 1.0 +
        2 * AA_Aa * 1.0 +
        2 * AA_aa * 1.0 +
        2 * Aa_Aa * 0.75 +
        2 * Aa_aa * 0.5 +
        2 * aa_aa * 0.0
    )

    return expected_dominant

# Example usage:
# Input: 1 2 3 4 5 6
# Output: 26


if __name__ == "__main__":
    input_data = input("Enter the six integers separated by spaces: ").strip()
    genotypes = list(map(int, input_data.split()))
    if len(genotypes) != 6:
        raise ValueError("Please provide exactly six integers.")
    result = expected_offsprings(genotypes)
    print(f"Expected number of dominant phenotype offspring: {result}")
