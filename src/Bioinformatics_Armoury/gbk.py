# URL: https://rosalind.info/problems/gbk/
# GenBank Introduction

# Given:
# Return:

from Bio import Entrez


def count_genbank_entries(genus, start_date, end_date):
    # provide your email when accessing NCBI API
    Entrez.email = "your email id"  # replace with your email id

    query = f'"{genus}"[Organism] AND ("{start_date}"[PDAT] : "{end_date}"[PDAT])'

    with Entrez.esearch(db="nucleotide", term=query) as handle:
        record = Entrez.read(handle)

    return int(record["Count"])


if __name__ == "__main__":
    genus = input("Please enter the Genus: ").strip()
    start_date = input("Provide the start date (YYYY/MM/DD): ").strip()
    end_date = input("Provide the end date (YYYY/MM/DD): ").strip()

    count = count_genbank_entries(genus, start_date, end_date)
    print(count)
