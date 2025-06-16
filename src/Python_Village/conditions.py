# Conditions & Loops
# URL: https://rosalind.info/problems/ini4/
# Given: Two positive integers a and b (a < b < 10000).
# Return: The sum of all odd integers from a through b, inclusively.

def sum_of_odds(a, b):
    # Initialize the sum
    total_sum = 0

    # Iterate through the range from a to b
    for i in range(a, b + 1):
        # Check if the number is odd
        if i % 2 != 0:
            total_sum += i

    return total_sum


if __name__ == "__main__":
    input_path = input("Enter the input file path: ")
    with open(input_path, 'r') as file:
        # Read the integers a and b from the file
        a, b = map(int, file.readline().strip().split())

    result = sum_of_odds(a, b)
    print(result)
