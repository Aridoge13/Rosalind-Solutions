# Strings and Lists
# URL: https://rosalind.info/problems/ini3/
# Given: A string s of length at most 200 letters and four integers a, b, c, and d.
# Return: The slice of this string from indices a through b and c through d (with space in between), inclusively. In other words, we should include elements s[b] and s[d] in our slice.

def string_slice(s, a, b, c, d):
    # Extract the slices from the string
    slice1 = s[a:b+1]  # +1 because the end index is inclusive
    slice2 = s[c:d+1]  # +1 because the end index is inclusive

    # Return the two slices separated by a space
    return f"{slice1} {slice2}"


if __name__ == "__main__":
    input_path = input("Enter the input file path: ")
    with open(input_path, 'r') as file:
        # Read the first line for the string
        s = file.readline().strip()
        # Read the second line for the integers
        a, b, c, d = map(int, file.readline().strip().split())
    result = string_slice(s, a, b, c, d)
    print(result)
