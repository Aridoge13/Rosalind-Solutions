# Counting Point Mutations
# Problem Title: Point Mutation
# Problem ID: 1A
# URL: http://rosalind.info/problems/hamm/

def count_point_mutations(s1, s2):
    """
    Count the number of point mutations (differences) between two strings of equal length.

    Parameters:
    s1 (str): First string.
    s2 (str): Second string.

    Returns:
    int: Number of point mutations.
    """
    if len(s1) != len(s2):
        raise ValueError("Strings must be of equal length.")

    return sum(1 for a, b in zip(s1, s2) if a != b)

# Example Sequence:
# s1 = "GAGCCTACTAACGGGAT"
# s2 = "CATCGTAATGACGGCCT"


if __name__ == "__main__":
    input_path = input("Enter the path to the input file: ")
    with open(input_path, 'r') as file:
        s1 = file.readline().strip()
        s2 = file.readline().strip()
    result = count_point_mutations(s1, s2)
    print(result)
