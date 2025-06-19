# Completing a Tree
# URL: https://rosalind.info/problems/tree/

# Given: A positive integer n (n≤1000) and an adjacency list corresponding to a graph on n nodes that contains no cycles.
# Return: The minimum number of edges that can be added to the graph to produce a tree.

def completing_tree(n, edges):
    # A tree with n nodes has exactly n-1 edges
    # The number of edges needed to complete the tree is (n - 1) - current_edges
    current_edges = len(edges)
    return (n - 1) - current_edges


def read_input(file_path):
    """
    Reads the input file and returns the number of nodes and the edges.
    """
    with open(file_path, 'r') as f:
        n = int(f.readline().strip())
        edges = []
        for line in f:
            if line.strip():
                edge = list(map(int, line.strip().split()))
                edges.append(edge)
    return n, edges


if __name__ == "__main__":
    input_path = input("Enter the path to the input file: ")
    n, edges = read_input(input_path)
    result = completing_tree(n, edges)
    print(result)
