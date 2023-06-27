 import Pkg
Pkg.add("Combinatorics")
using Combinatorics
using LinearAlgebra
include("helpers.jl")
include("greedy_max_cut.jl")

"""
    random_y(x,t,n)
    x: normalized min eigenvector
    t: randomly chosen value to partition vertices
    n: length of x

    returns a vector p indicating how t partitions the vertices in the cut

"""
function random_y(x,t,n)
    p = zeros(n)
    for i = 1:n
        if x[i] <= -sqrt(t)
            p[i] = -1
        elseif x[i] >= sqrt(t)
            p[i] = 1
        end
    end
    return p
end


"""
    trevisan_max_cut(adj_matrix, adj_list, active_verts, n, numt)
    adj_matrix : the adjacency matrix of the input graph
    adj_list : the adjacency lists of the input graph
    active_verts : vertices that have not been partitioned in a previous iteration
    n : the number of vertices
    numt: the number of times to run the max cut algorithm

returns, for a given graph in adjacency matrix form `adj_matrix`,
the 0.531-approximation max cut value and its corresponding vertex
set partition.

This is Recursive-Spectral-Cut from the Trevisan paper
"""
function trevisan_max_cut(adj_matrix, adj_list, active_verts, n, numt)

    D = adj_matrix_to_degree_matrix(adj_matrix,n)
    A = adj_matrix
    neg_sqrt_D = zeros(n,n)
    [neg_sqrt_D[i,i] = D[i,i]==0 ? 0 : D[i,i]^(-1/2) for i ∈ 1:n]
    x = real.(smallest_eigenvector(neg_sqrt_D * A * neg_sqrt_D)) # sometimes theres a 0.0im part
    max_x = maximum(abs.(x))

    for i ∈ 1:n
        x[i] = x[i]/max_x
    end

    current_max = 0
    M = 0
    C = 0
    X = 0
    final_y = zeros(n)

    #change this line to test more than 1 t value per iteration
    for k = 1:1

        t = rand(MersenneTwister())
        y = random_y(x,t,n)
        m1 = 0
        m2 = 0
        c = 0
        xx = 0
        for i ∈ 1:n
            for j ∈ adj_list[i]
                if j > i
                    m1 += adj_matrix[i,j]
                    if y[i] == 0 && y[j] == 0
                        m2 += adj_matrix[i,j]
                    elseif y[i]*y[j] == -1
                        c += adj_matrix[i,j]
                    elseif abs(y[i] + y[j]) == 1
                        xx += adj_matrix[i,j]
                    end
                end
            end
        end

        m = (m1-m2)

        if c + xx/2 - m/2 > current_max
            current_max = c + xx/2 - m/2
            M = m
            C = c
            X = xx
            final_y = y
        end
    end

    y = final_y


    if (C + X/2 <= M/2)
        greedy_sol = greedy_max_cut(adj_matrix, active_verts, adj_list, n)
        cut = greedy_sol[1]
        return (greedy_sol[2],cut)
    else
        L = filter(i -> y[i] == -1 && i ∈ active_verts, 1:n)
        R = filter(i -> y[i] == 1 && i ∈ active_verts, 1:n)
        V_prime = filter(i -> y[i] == 0 && i ∈ active_verts, 1:n)

        if !isempty(V_prime)
            (cut, cut_val) = trevisan_max_cut(induced_subgraph(adj_matrix,V_prime),adj_list, V_prime, n, numt)
        else
            cut = []
        end

        cut_val_1 = cut_value(adj_matrix, union(cut,L), adj_list)
        cut_val_2 = cut_value(adj_matrix, union(cut,R), adj_list)

        if cut_val_1 > cut_val_2
            return union(cut,L), cut_val_1
        else
            return union(cut,R), cut_val_2
        end
    end
end
