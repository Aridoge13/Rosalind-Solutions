# Intro to Python Dictionaries
# URL: https://rosalind.info/problems/ini6
# Given: A string s of length at most 10000 letters.
# Return: The number of occurrences of each word in s, where words are separated by spaces. Words are case-sensitive, and the lines in the output can be in any order.

# Sample Input: We tried list and we tried dicts also we tried Zen
# Sample Output:
# and 1
# We 1
# tried 3
# dicts 1
# list 1
# we 2
# also 1
# Zen 1


def count_word_occurrences(s):
    # Split the string into words
    words = s.split()

    # Create a dictionary to count occurrences
    word_count = {}

    # Count occurrences of each word
    for word in words:
        if word in word_count:
            word_count[word] += 1
        else:
            word_count[word] = 1

    return word_count


def print_word_counts(word_count):
    # Print the word counts
    for word, count in word_count.items():
        print(f"{word} {count}")


if __name__ == "__main__":
    input_path = input("Enter the input file path: ")
    with open(input_path, 'r') as file:
        s = file.read().strip()  # Read the entire file and strip any leading/trailing whitespace
    word_count = count_word_occurrences(s)
    print_word_counts(word_count)
