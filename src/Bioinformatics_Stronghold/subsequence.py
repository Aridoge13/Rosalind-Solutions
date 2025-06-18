# Longest Increasing Subsequence
# URL: https://rosalind.info/problems/lgis/


def longest_subsequence(perm, increasing=True):
    n = len(perm)
    dp = [1] * n
    prev = [-1] * n

    for i in range(n):
        for j in range(i):
            if (increasing and perm[i] > perm[j]) or (not increasing and perm[i] < perm[j]):
                if dp[j] + 1 > dp[i]:
                    dp[i] = dp[j] + 1
                    prev[i] = j

    # Find index of max length
    max_len = max(dp)
    index = dp.index(max_len)

    # Backtrack
    seq = []
    while index != -1:
        seq.append(perm[index])
        index = prev[index]

    return seq[::-1]  # reverse the sequence


if __name__ == "__main__":
    input_path = input("Enter the input file path: ")
    with open(input_path, 'r') as input_file:
        lines = input_file.read().strip().splitlines()
    n = int(lines[0])
    perm = list(map(int, lines[1].split()))

    lis_result = longest_subsequence(perm, increasing=True)
    lds_result = longest_subsequence(perm, increasing=False)

    output_path = input("Enter the output file path: ")
    with open(output_path, 'w') as output_file:
        output_file.write(' '.join(map(str, lis_result)) + '\n')
        output_file.write(' '.join(map(str, lds_result)) + '\n')
