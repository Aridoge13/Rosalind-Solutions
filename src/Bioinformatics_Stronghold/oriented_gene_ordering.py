# Enumerating Oriented Gene Orderings
# URL: https://rosalind.info/problems/sign/


# Given: A positive integer n (n ≤ 6).
# Return: The total number of signed permutations of length n, followed by a list of all such permutations(you may list the signed permutations in any order).


from itertools import permutations


def signed_permutations(n):
    # Generate all permutations of the numbers from 1 to n
    perms = permutations(range(1, n + 1))

    # Create a list to hold the signed permutations
    signed_perms = []

    # For each permutation, create all combinations of signs
    for perm in perms:
        for signs in range(2 ** n):
            signed_perm = []
            for i in range(n):
                if (signs >> i) & 1:
                    signed_perm.append(perm[i])  # positive
                else:
                    signed_perm.append(-perm[i])  # negative
            signed_perms.append(tuple(signed_perm))

    return len(signed_perms), signed_perms


if __name__ == "__main__":
    n = int(input("Enter a positive integer n (n ≤ 6): "))
    count, permutations_list = signed_permutations(n)
    print(count)
    for perm in permutations_list:
        print(" ".join(map(str, perm)))
