using JuMP
import LinearAlgebra
import Random
import SCS
import Test
using TimerOutputs
using BenchmarkTools

"""
best_hyperplane(cut,n,V,weights,adj_list)
cut : array to fill with entries indicating which side of the cut each vertex is in
n : number of vertices
V : solution from SDP solver
weights: edge weights
adj_list : adjacency list of input graph

function that tests 100 random hyperplanes to partition vertices based on the solution from the SDP
    returns the best of these cuts
"""
function best_hyperplane(cut,n,V,weights,adj_list)

    M = MersenneTwister()
    rand(M)
    r = rand(M, 100, n)

    max_cut_value = 0
    max_cut = Array{Int}(undef,n)
    for j in 1:100
        A = []

        for i in 1:n
            if LinearAlgebra.dot(r[j,:], V[:, i]) <= 0
                cut[i] = -1
                push!(A,i)
            else
                cut[i] = 1
            end
        end

        #compute edge by edge whether cut value of 2 endpoints is different and add weight to cut
        temp_cut_val = cut_value(weights, A, adj_list)


        if  temp_cut_val > max_cut_value
            max_cut_value = temp_cut_val
            max_cut = cut
        end
    end
    return max_cut_value, max_cut
end

"""
solve(num_vertex, weights)
num_vertex : number of vertices
weights : input edge weights

returns a matrix V resulting from solving the Max-Cut SDP
note: the cut is not approximated in this step/function
"""
function solve(num_vertex, weights)
    # Calculate the (weighted) Lapacian of the graph: L = D - W.
    laplacian = LinearAlgebra.diagm(0 => weights * ones(num_vertex)) - weights
    # Solve the SDP relaxation
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, X[1:num_vertex, 1:num_vertex], PSD)
    @objective(model, Max, 1 / 4 * LinearAlgebra.dot(laplacian, X))
    @constraint(model, LinearAlgebra.diag(X) .== 1)
    optimize!(model)
    # Compute the Cholesky factorization of X, i.e., X = V^T V.
    opt_X = LinearAlgebra.Hermitian(value.(X), :U)  # Tell Julia its PSD.
    factorization = LinearAlgebra.cholesky(opt_X, Val(true); check = false)
    V = (factorization.P * factorization.L)'
    return V
end

"""
solve_max_cut_sdp(num_vertex,weights, adj_list)
num_vertex : number of vertices
weights : input edge weights
adj_list : adjacency list of input graph

returns the best cut and cut value using the SDP approximation algorithm and 100 random hyperplanes
"""
function solve_max_cut_sdp(num_vertex, weights, adj_list)

    V = solve(num_vertex, weights)
    # Normalize columns.
    for i in 1:num_vertex
        V[:, i] ./= LinearAlgebra.norm(V[:, i])
    end
    # Generate random vector on unit sphere.

    cut = Array{Int}(undef,num_vertex)
    max_cut_results= best_hyperplane(cut,num_vertex, V, weights, adj_list)

    return max_cut_results[1], max_cut_results[2]
end
