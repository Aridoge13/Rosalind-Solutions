# Transcribing DNA to RNA

def transcribe_dna_to_rna(dna_sequence):
    """
    Transcribes a DNA sequence to an RNA sequence by replacing 'T' with 'U'.

    Parameters:
    dna_sequence (str): The DNA sequence to be transcribed.

    Returns:
    str: The transcribed RNA sequence.
    """
    return dna_sequence.replace('T', 'U')


if __name__ == "__main__":
    input_path = input("Enter the path to your input file here: ")
    output_path = input("Enter the path to your output file here: ")
    with open(input_path, 'r') as file:
        dna_sequence = file.read().strip()
    rna_sequence = transcribe_dna_to_rna(dna_sequence)
    with open(output_path, 'w') as file:
        file.write(rna_sequence)
