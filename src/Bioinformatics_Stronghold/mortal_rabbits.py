# Mortal Fibonacci Rabbits
# Wabbit Season
# URL: https://rosalind.info/problems/fibd/
# Given: Positive integers n and m (n ≤ 100, m ≤ 20).
# Return: The total number of pairs of rabbits that will remain after the n-th month, if all rabbits live for m months.

def mortal_fibonacci_rabbits(n, m):
    # Each index represents the number of rabbit pairs of a given age
    ages = [0] * m
    ages[0] = 1  # Start with 1 pair of newborns

    for month in range(1, n):
        # Rabbits that are at least 1 month old can reproduce
        new_born = sum(ages[1:])
        # Age the population: shift right, drop oldest, add new_born at ages[0]
        ages = [new_born] + ages[:-1]
    return sum(ages)


if __name__ == "__main__":
    n, m = map(int, input("Enter n and m: ").split())
    result = mortal_fibonacci_rabbits(n, m)
    # Total pairs of rabbits will be printed to the terminal.
    print("Total pairs of rabbits:", result)
