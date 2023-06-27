include("helpers.jl")
using LinearAlgebra

"""
split_cuts(A,n,adj_list)
A : adjacency matrix of input graph
n : number of vertices
adj_list : adjacency list of input graph

returns the best value of all n-1 possible split cuts (sort the min eigenvector and consider partitions given by splitting on either side of a chosen index)
"""
function split_cuts(A, n,adj_list)
    D = adj_matrix_to_degree_matrix(A,n)
    neg_sqrt_D = zeros(n,n)
    [neg_sqrt_D[i,i] = D[i,i]==0 ? 0 : (D[i,i]+0im)^(-1/2) for i ∈ 1:n]
    x = real.(smallest_eigenvector(neg_sqrt_D * A * neg_sqrt_D))

    p = sortperm(x)

    max_val = 0
    max_cut = []
    for j = 1:n-1
        cut, val = local_search(A, adj_list, p[1:j], n, cut_value(A, p[1:j], adj_list))
        if val > max_val
            max_val = val
            max_cut = p[1:j]
        end
    end
    return max_val

end
