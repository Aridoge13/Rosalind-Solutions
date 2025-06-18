# Calculating Protein Mass
# URL: https://rosalind.info/problems/prtm/

# Given: A Protein string P of length at most 1000 amino acids.
# Return: The total weight of P. Consult the monoisotopic mass table.

from collections import defaultdict


def calculate_protein_mass(protein):
    # Monoisotopic mass table for amino acids
    mass_table = {
        'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259,
        'F': 147.06841, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406,
        'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
        'P': 97.05276, 'Q': 128.05858, 'R': 156.10111, 'S': 87.03203,
        'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333
    }

    # Calculate the total mass of the protein string
    total_mass = sum(mass_table[aa] for aa in protein)

    return total_mass


if __name__ == "__main__":
    protein = input(
        "Enter the primary structure of your protein here:").strip()
    total_mass = calculate_protein_mass(protein)
    print(total_mass)

# Example usage:
# If the input is "SKADYEK", the output will be 821.392
# Explanation: The total mass is calculated as follows:
# S (87.03203) + K (128.09496) + A (71.03711) + D (115.02694) + Y (163.06333) + E (129.04259) + K (128.09496) = 821.392
