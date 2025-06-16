# Wascally Wabbits
# Problem Title: Rabbits and Recurrence Relations
# URL: http://rosalind.info/problems/fib/
def rabbits(n, k):
    """
    Calculate the number of rabbit pairs after n months, given that each pair produces k new pairs every month.

    :param n: Total number of months
    :param k: Number of new pairs produced by each pair every month
    :return: Total number of rabbit pairs after n months
    """
    if n == 1 or n == 2:
        return 1

    prev = 1  # F(1)
    curr = 1  # F(2)

    for month in range(3, n + 1):
        next_pair = curr + k * prev
        prev = curr
        curr = next_pair

    return curr


if __name__ == "__main__":
    input_data = input("Enter n and k (space-separated): ")
    n, k = map(int, input_data.split())
    result = rabbits(n, k)
    print(result)
# Example usage:
# n, k = 5, 3
# print(rabbits(n, k))  # Output: 19
