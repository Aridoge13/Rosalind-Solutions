# Finding a shared motif
# URL: https://rosalind.info/problems/lcsm/
# Given: A collection of k(<=100) DNA strings of length at most 1 kbp each in FASTA format.
# A longest common substring of the collection. (If multiple solutions exist, you may return any single solution.)

def shared_motif(fasta_strings):
    """
    Find the longest common substring in a collection of DNA strings.

    :param fasta_strings: List of DNA strings.
    :return: The longest common substring.
    """
    if not fasta_strings:
        return ""
    # Start with the first string as the initial candidate
    shortest_string = min(fasta_strings, key=len)
    longest_motif = ""
    # Check all substrings of the shortest string
    for length in range(len(shortest_string), 0, -1):
        for start in range(len(shortest_string) - length + 1):
            candidate = shortest_string[start:start + length]
            if all(candidate in s for s in fasta_strings):
                return candidate  # Return the first found longest motif

    return longest_motif  # If no common motif is found


if __name__ == "__main__":
    input_path = input("Enter the path to the FASTA file: ").strip()
    with open(input_path, 'r') as file:
        input_data = file.readlines()
    fasta_strings = []
    current_string = ""
    for line in input_data:
        if line.startswith(">"):
            if current_string:
                fasta_strings.append(current_string)
                current_string = ""
        else:
            current_string += line.strip()
    if current_string:
        fasta_strings.append(current_string)
    longest_motif = shared_motif(fasta_strings)
    print(f"Longest common substring: {longest_motif}")
