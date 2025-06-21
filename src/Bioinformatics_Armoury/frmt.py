# Data Formats
# URL: https://rosalind.info/problems/frmt/

# Given: A collection of n(n≤10) GenBank entry IDs.
# Return: The shortest of the strings associated with the IDs in FASTA format.


from Bio import SeqIO
from Bio import Entrez


def fetch_shortest_fasta(ids):
    Entrez.email = "aritra.mukherjee98@gmail.com"  # Use your actual email
    with Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text") as handle:
        records = list(SeqIO.parse(handle, "fasta"))

    # Find the shortest sequence
    shortest_record = min(records, key=lambda record: len(record.seq))
    return shortest_record


if __name__ == "__main__":
    # Accept input from user or file
    ids = input("Enter GenBank IDs separated by space: ").strip().split()

    shortest = fetch_shortest_fasta(ids)

    # Print FASTA format
    print(f">{shortest.description}")
    print(shortest.seq)
