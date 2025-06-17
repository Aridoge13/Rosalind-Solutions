# Overlap Graphs
# URL: http://rosalind.info/problems/grph/
# Given: A collection of DNA strings in FASTA format having total length at most 10 kbp.
# Return: The adjacency list corresponding to O3. You may return edges in any order.

from Bio import SeqIO


def overlap_graphs(fasta_file, k=3):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    seqs = [(record.id, str(record.seq)) for record in records]
    edges = []
    for id1, seq1 in seqs:
        for id2, seq2 in seqs:
            if id1 != id2 and seq1[-k:] == seq2[:k]:
                edges.append((id1, id2))
    return edges


if __name__ == "__main__":
    fasta_path = input("Enter the path to the FASTA file: ").strip()
    edges = overlap_graphs(fasta_path, k=3)
    output_path = input("Enter the path to save the output: ").strip()
    with open(output_path, "w") as output_file:
        for id1, id2 in edges:
            output_file.write(f"{id1} {id2}\n")
