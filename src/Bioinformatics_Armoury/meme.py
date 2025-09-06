# New Motif Discovery
# URL: https://rosalind.info/problems/meme/


# Given: A set of protein strings in FASTA format that share some motif with minimum length 20.
# Return: Regular expression for the best-scoring motif.

from Bio import SeqIO


def read_fasta_files(fasta_file):
    return [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]


# Sample Input:
# >Rosalind_7142
# PFTADSMDTSNMAQCRVEDLWWCWIPVHKNPHSFLKTWSPAAGHRGWQFDHNFFVYMMGQ
# FYMTKYNHGYAPARRKRFMCQTFFILTFMHFCFRRAHSMVEWCPLTTVSQFDCTPCAIFE
# WGFMMEFPCFRKQMHHQSYPPQNGLMNFNMTISWYQMKRQHICHMWAEVGILPVPMPFNM
# SYQIWEKGMSMGCENNQKDNEVMIMCWTSDIKKDGPEIWWMYNLPHYLTATRIGLRLALY
# >Rosalind_4494
# VPHRVNREGFPVLDNTFHEQEHWWKEMHVYLDALCHCPEYLDGEKVYFNLYKQQISCERY
# PIDHPSQEIGFGGKQHFTRTEFHTFKADWTWFWCEPTMQAQEIKIFDEQGTSKLRYWADF
# QRMCEVPSGGCVGFEDSQYYENQWQREEYQCGRIKSFNKQYEHDLWWCWIPVHKKPHSFL
# KTWSPAAGHRGWQFDHNFFSTKCSCIMSNCCQPPQQCGQYLTSVCWCCPEYEYVTKREEM
# >Rosalind_3636
# ETCYVSQLAYCRGPLLMNDGGYGPLLMNDGGYTISWYQAEEAFPLRWIFMMFWIDGHSCF
# NKESPMLVTQHALRGNFWDMDTCFMPNTLNQLPVRIVEFAKELIKKEFCMNWICAPDPMA
# GNSQFIHCKNCFHNCFRQVGMDLWWCWIPVHKNPHSFLKTWSPAAGHRGWQFDHNFFQMM
# GHQDWGTQTFSCMHWVGWMGWVDCNYDARAHPEFYTIREYADITWYSDTSSNFRGRIGQN

def find_shared_motif(sequences, min_length=20):
    if not sequences:
        return ""
    shortest_seq = min(sequences, key=len)
    max_len = len(shortest_seq)
    for length in range(max_len, min_length - 1, -1):
        for start in range(0, max_len - length + 1):
            motif = shortest_seq[start:start + length]
            if all(motif in seq for seq in sequences):
                return motif
    return ""


def motif_to_regex(motif):
    return motif

# Sample output: DLWWCWIPVHK[NK]PHSFLKTWSPAAGHRGWQFDHNFF


def consensus_motif(sequences):
    from itertools import zip_longest
    motif = ""
    for chars in zip(*sequences):
        unique = set(chars)
        if len(unique) == 1:
            motif += unique.pop()
        else:
            motif += "[" + "".join(sorted(unique)) + "]"
    return motif


if __name__ == "__main__":
    fasta_file = input("Enter path to FASTA file: ").strip()
    sequences = read_fasta_files(fasta_file)
    print("Number of sequences:", len(sequences))
    print("Lengths:", [len(seq) for seq in sequences])
    print("Shortest sequence length:", len(
        min(sequences, key=len)) if sequences else 0)
    motif = consensus_motif(sequences)
    if motif:
        regex = motif_to_regex(motif)
        print(f"Best-scoring motif (regex): {regex}")
    else:
        print("No shared motif of minimum length 20 found.")
