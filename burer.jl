using JuMP, Ipopt, Calculus, Random
include("helpers.jl")
include("linesearch.jl")
include("bmz_linesearch.jl")

function f(theta)
    n = length(theta)
    return (1/2)*sum(sum(A[i,j]*cos(theta[i]-theta[j]) for j in 1:n) for i in 1:n)
end

function grad_f(theta)
    n = length(theta)
    return [sum(A[k,j]*sin(theta[k] - theta[j]) for k in 1:n) for j in 1:n]
end

function hess_f(theta)
    n = length(theta)
    hess = zeros(n,n)
    for i = 1:n
        for j = 1:n
            if i == j
                hess[i,i] = -sum(A[k,i]*cos(theta[k]-theta[i]) for k in 1:n) + A[i,i]*cos(theta[i] - theta[i])
            else
                hess[i,j] = A[i,j]*cos(theta[i] - theta[j])
            end
        end
    end
    return hess
end

A = []
adj_list = []


function procedure_cut(theta)
    n = length(theta)
    best_cut = []
    alpha = 0
    gamma = -Inf
    push!(theta, 2*π)
    for i=1:n+1
        cut = generate_cut(theta[i], theta)
        cut_val = cut_value(A, cut, adj_list)
        if cut_val > gamma
            gamma = cut_val
            best_cut = cut
        end
    end
    return best_cut, gamma
end

function generate_cut(alpha, theta)
    n = length(theta) - 1
    if alpha <= π
        cut = filter(i -> theta[i] < alpha + π && theta[i] >= alpha, 1:n)
    else
        cut = filter(i -> theta[i] < alpha  && theta[i] >= alpha - π, 1:n)
    end
    return cut
end

function windmill_cut(alpha, k, theta)
    n = length(theta) -1
    cut = []
    for j = 1:k
        temp_cut = filter(i -> theta[i] >= alpha + (2*j-2)*(π/k) && theta[i] < alpha + (2*j-1)*(π/k), 1:n)
        cut = union(cut, temp_cut)
    end
    cut_val = cut_value(A, cut, adj_list)
    return cut, cut_val
end

function swap_one(cut,n, val, change, tolerance)
    swap = true
    while swap
        swap = false
        for i =1:n
            if change[i] > tolerance
                val = val + change[i]
                cut,change = updatecut(i,cut,change)
                swap = true
            end
        end
    end
    return cut, val
end

function updatecut(i,cut,change)
    if i ∈ cut
        for j ∈ adj_list[i]
            if j ∈ cut
                change[j] = change[j] - 2*A[i,j]
            else
                change[j] = change[j] + 2*A[i,j]
            end
        end
        cut = filter(p -> p != i, cut)
    else
        for j ∈ adj_list[i]
            if j ∉ cut
                change[j] = change[j] - 2*A[i,j]
            else
                change[j] = change[j] + 2*A[i,j]
            end
        end
        cut = union(cut, [i])
    end
    return cut, change
end


function burer_heuristic(adj_matrix, n, adj, N)

    global A = adj_matrix
    global adj_list = adj
    k= 0
    gamma = -Inf
    best_cut = []


    temp = rand(n)
    theta_init = 2*π.*temp

    w1norm = sum(sum(A[i,j] for j in 1:n) for i in 1:n)

    p = π/n
    while k <= N

    opt_theta = linesearch(theta_init, A, w1norm, n)

        for i =1:n
            if opt_theta[i] > (2*π)
                opt_theta[i] = opt_theta[i] % (2*π)
            elseif opt_theta[i] < 0
                while opt_theta[i] < 0
                    opt_theta[i] = opt_theta[i] + 2*π
                end
            end
        end
        cut, cut_value = procedure_cut(opt_theta)
        cut, cut_value = local_search(adj_matrix, adj, cut, n, cut_value)

        if cut_value > gamma
            gamma = cut_value
            best_cut = cut
            k = 0
        else
            k = k + 1
        end
        y = zeros(n)
        for i = 1:n
            if i ∈ cut
                y[i] = π
            end
        end

        theta_init = [y[i] + rand([-π,π])*.2 for i in 1:n]
    end

    return gamma, best_cut

end
