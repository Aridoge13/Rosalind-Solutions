# Computing GC Content
# Rosalind ID: GC
# URL: http://rosalind.info/problems/gc/
from Bio import SeqIO


def compute_gc_content(fasta_file):
    max_gc_content = 0.0
    max_gc_id = ""

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        gc_content = (seq.count('G') + seq.count('C')) / len(seq) * 100

        if gc_content > max_gc_content:
            max_gc_content = gc_content
            max_gc_id = record.id
    return max_gc_id, max_gc_content


if __name__ == "__main__":
    fasta_file = input("Enter the path to the FASTA file: ").strip()
    # Read the FASTA file and compute the GC content
    # and the ID of the sequence with the highest GC content
    max_id, max_gc = compute_gc_content(fasta_file)
    print(f"{max_id}\n{max_gc:.6f}")
