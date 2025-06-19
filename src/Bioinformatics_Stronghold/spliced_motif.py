# Finding a spliced motif
# URL: https://rosalind.info/problems/sseq/

# Given: Two DNA strings s and t (each of length at most 1 kbp) in FASTA format.
# Return: One collection of indices of s in which the symbols of t appear as a subsequence of s. If multiple solutions exist, you may return any one.

def parse_fasta(file_path):
    with open(file_path, 'r') as f:
        lines = f.read().strip().split('>')
        records = []
        for entry in lines:
            if not entry:
                continue
            parts = entry.split('\n')
            seq = ''.join(parts[1:]).replace(
                '\r', '').replace('\n', '').strip()
            records.append(seq.upper())
        return records


def find_spliced_motif(s, t):
    indices = []
    start = 0
    for char in t:
        pos = s.find(char, start)
        if pos == -1:
            raise ValueError(
                f"Character '{char}' not found after position {start}")
        indices.append(pos + 1)  # 1-based index
        start = pos + 1
    return indices


if __name__ == "__main__":
    fasta_path = input("Enter FASTA file path: ").strip()
    s, t = parse_fasta(fasta_path)

    print("DEBUG: s =", s)
    print("DEBUG: t =", t)

    result = find_spliced_motif(s, t)
    print(' '.join(map(str, result)))
