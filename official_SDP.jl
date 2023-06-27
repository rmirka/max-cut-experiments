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
function best_hyperplane(cut,n,V,weights,adj_list, num_vectors, search_each)

    best_A = []
    M = MersenneTwister()
    rand(M)
    r = rand(M, num_vectors, n)

    max_cut_value = 0
    max_cut = Array{Int}(undef,n)
    for j in 1:num_vectors
        A = []

        for i in 1:n
            if LinearAlgebra.dot(r[j,:], V[:, i]) <= 0
                cut[i] = -1
                push!(A,i)
            else
                cut[i] = 1
            end
        end

        #local_search every random vector
        if search_each
            A, temp_cut_val = local_search(weights,adj_list,A,n,cut_value(weights, A, adj_list))
        else
            temp_cut_val = cut_value(weights, A, adj_list)
        end

        if  temp_cut_val > max_cut_value
            max_cut_value = temp_cut_val
            max_cut = cut
            best_A = A
        end
    end

    #just local search current best vector
    if !search_each
        max_cut, max_cut_value = local_search(weights, adj_list, best_A, n, max_cut_value)
    end

    return max_cut_value, max_cut
end



#continue to try >= 1 new random vectors until time limit reached
function best_hyperplane_timed(cut,n,V,weights,adj_list,limit)

    t = time()

    best_A = []
    M = MersenneTwister()

    max_cut_value = 0
    max_cut = Array{Int}(undef,n)
    c = 0
    while c==0 || time() < t + limit
        c += 1
        r = rand(M,n)
        A = []

        for i in 1:n
            if LinearAlgebra.dot(r, V[:, i]) <= 0
                cut[i] = -1
                push!(A,i)
            else
                cut[i] = 1
            end
        end

        #local_search every random vector
        A, temp_cut_val = local_search(weights,adj_list,A,n,cut_value(weights, A, adj_list))


        if  temp_cut_val > max_cut_value
            max_cut_value = temp_cut_val
            max_cut = cut
            best_A = A
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

returns the best cut and cut value using the SDP approximation algorithm and 10 random hyperplanes with local search
"""
function solve_max_cut_sdp(num_vertex, weights, adj_list, time_limit, limit, num_vectors, search_each)

    t = time()
    V = solve(num_vertex, weights)
    # Normalize columns.
    for i in 1:num_vertex
        V[:, i] ./= LinearAlgebra.norm(V[:, i])
    end

    cut = Array{Int}(undef,num_vertex)
    if !time_limit
        max_cut_results= best_hyperplane(cut,num_vertex, V, weights, adj_list, num_vectors, search_each)
    else
        max_cut_results= best_hyperplane_timed(cut, num_vertex, V, weights, adj_list, limit-(time()-t))
    end

    return max_cut_results[1], max_cut_results[2]
end
