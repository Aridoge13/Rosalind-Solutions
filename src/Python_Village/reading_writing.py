# Working with Files
# Reading and Writing
# Given: A file containing at most 1000 lines.
# Return: A file containing all the even-numbered lines from the original file. Assume 1-based numbering of lines.
# URL: https://rosalind.info/problems/ini5/

def extract_even_lines(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line_number, line in enumerate(infile, start=1):
            # Check if the line number is even
            if line_number % 2 == 0:
                outfile.write(line)


if __name__ == "__main__":
    input_path = input("Enter the input file path: ")
    output_path = input("Enter the output file path: ")
    extract_even_lines(input_path, output_path)
    print(f"Even lines extracted to {output_path}")
