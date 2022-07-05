import Pkg
Pkg.add("Combinatorics")
using Combinatorics
include("helpers.jl")


"""
    greedy_max_cut(adj_matrix)

returns, for a given graph in adjacency matrix form `adj_matrix`,
the 0.5-approximation max cut value and its corresponding vertex
set partition.

This is a greedy algorithm and runs in time linear to the number
of vertices.
"""
function greedy_max_cut(adj_matrix, active_verts, adj_list, n)
    A = []
    B = []
    value = 0
    for vertex in active_verts
        
        relevant_row = adj_matrix[vertex, :]
        A_neighbors_weight = sum([relevant_row[i] for i ∈ filter(i -> relevant_row[i] != 0 && i in A, adj_list[vertex])])
        B_neighbors_weight = sum([relevant_row[i] for i ∈ filter(i -> relevant_row[i] != 0 && i in B, adj_list[vertex])])

        if A_neighbors_weight > B_neighbors_weight
            B = union(B, [vertex])
            value += A_neighbors_weight
        else
            A = union(A, [vertex])
            value += B_neighbors_weight
        end

    end

    return (value, A)
end
