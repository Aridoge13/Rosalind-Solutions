# Translating RNA to Protein
# URL: http://rosalind.info/problems/rna/
# Translating RNA to Protein
from Bio.Seq import Seq


def translate_rna_to_protein(rna_sequence):
    """
    Translate an RNA sequence into a protein sequence.
    :param rna_sequence: A string representing the RNA sequence.
    :return: A string representing the translated protein sequence.
    """
    rna_seq = Seq(rna_sequence)
    # Translate RNA to protein using BioPython's Seq module
    protein_seq = rna_seq.translate()
    return str(protein_seq)


# Provide the path to the input file containing the RNA sequence
if __name__ == "__main__":
    input_path = input("Enter the path to the input file: ")
    with open(input_path, 'r') as file:
        rna_sequence = file.readline().strip()
    result = translate_rna_to_protein(rna_sequence)
    print(result)

# Example Input:
# AUGGCCUAAUGG
# Example Output:
# MAAW
