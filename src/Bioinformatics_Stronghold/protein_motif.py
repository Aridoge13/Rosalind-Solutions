# Finding a Protein Motif
# URL: https://rosalind.info/problems/mprt/

# Given: At most 15 UniProt Protein Database access IDs.
# Return: For each protein possessing the N-glycosylation motif, output its given access ID followed by a list of locations in the protein string where the motif can be found.


import requests
import re


def find_protein_motif(protein_ids):
    motif = re.compile(r'N[^P][ST][^P]')
    results = []
    for protein_id in protein_ids:
        accession = protein_id.split('_')[0]
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
        response = requests.get(url)
        if response.status_code != 200:
            continue
        # Parse FASTA
        seq = ''.join(line.strip()
                      for line in response.text.splitlines() if not line.startswith('>'))
        # Find all motif matches (overlapping, 1-based)
        positions = [
            m.start() + 1 for m in re.finditer(r'(?=(N[^P][ST][^P]))', seq)]
        if positions:
            results.append((protein_id, positions))
    return results


# Example usage with sample protein IDs (replace with actual IDs or remove these lines if not needed)
# protein_ids = ['P07204_TRBM_HUMAN', 'P20840_SAG1_YEAST']
# results = find_protein_motif(protein_ids)
# for pid, positions in results:
#     print(pid, *positions)


if __name__ == "__main__":
    input_data = input("Enter the path to the file containing protein IDs: ")
    with open(input_data, 'r') as file:
        input_data = file.read().strip()
    protein_ids = input_data.split()
    output_path = input("Enter the output file path: ")
    results = find_protein_motif(protein_ids)
    with open(output_path, 'w') as file:
        for result in results:
            # result is a tuple: (protein_id, [positions])
            protein_id, positions = result
            line = protein_id + ' ' + ' '.join(map(str, positions))
            file.write(line + "\n")
