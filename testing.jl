include("trevisan_max_cut.jl")
include("official_SDP.jl")
include("helpers.jl")
include("split_cuts.jl")
using MatrixMarket
using LightGraphs
using Random
using SparseArrays
using TSPLIB

"""
simple_spectral(A,n,adj_list)
A : adjacency matrix of input graph
n : number of vertices
adj_list : adjacency list of input graph

returns the cut determined by splitting the minimum eigenvalue at 0
"""
function simple_spectral(A,n,adj_list)
    D = adj_matrix_to_degree_matrix(A,n)
    neg_sqrt_D = zeros(n,n)
    [neg_sqrt_D[i,i] = D[i,i]==0 ? 0 : (D[i,i]+0im)^(-1/2) for i ∈ 1:n]
    x = real.(smallest_eigenvector(neg_sqrt_D * A * neg_sqrt_D))
    L = filter(i -> x[i] <0, 1:n)
    return cut_value(A, L, adj_list)
end

"""
best_trevisan(A, adj, size, iter)
A : adjacency matrix of input graph
adj : adjacency list of input graph
size : number of vertices
iter : number of times to run trevisan's spectral algorithm

returns the best cut and value of 'iter' runs of trevisan's spectral algorithm
"""
function best_trevisan(A, adj, size, iter)
    max = 0
    max_cut = []
    for i =1:iter
        val, cut = trevisan_max_cut(A, adj, collect(1:size), size, size)
        if val > max
            max = val
            max_cut = cut
        end
    end
    return max, max_cut
end

#run these to prep the '@timed' function without having to process the function
@timed best_trevisan([0.0 1.0; 1.0 0.0], [[2],[1]], 2, 5)
@timed solve_max_cut_sdp(2, [0.0 1.0; 1.0 0.0], [[2],[1]])
@timed simple_spectral([0.0 1.0; 1.0 0.0],2, [[2],[1]])
@timed greedy_max_cut([0.0 1.0; 1.0 0.0], collect(1:2), [[2],[1]], 2)
@timed split_cuts([0.0 1.0; 1.0 0.0], 2, [[2][1]])


files = [:att48, :bayg29, :bays29, :berlin52, :bier127, :brazil58, :brg180, :ch130, :ch150, :d198, :eil101, :gr120, :gr137, :gr202, :gr96, :kroA100,:a280]

#random graphs with 50,100, 200, 350, 500 vertices and edge probs of .1, .25, .5, .75
sizes = [50,100,200,350,500]
p = [.1,.25,.5,.75]


spectral_times = []
spectral_cut_vals = []
trev_times = []
trev_cut_vals = []
sdp_times = []
sdp_cut_vals = []
greedy_times = []
greedy_cut_vals = []
split_times = []
split_cut_vals = []

#test average values for random graphs
# for i=1:4
#     for j = 1:4
#         adj = []
#         A = zeros(sizes[i],sizes[i])
#         G = erdos_renyi(sizes[i], p[j])
#         for k=1:sizes[i]
#             nk = neighbors(G,k)
#             push!(adj, nk)
#             for l in nk
#                 A[k,l] = 1
#             end
#         end
#
#         # for k = 1:6
#         #     avg_time = 0
#         #     avg_value = 0
#         #     for l = 1:10
#         #         t1 = @timed trevisan_max_cut(A, adj, collect(1:sizes[i]), sizes[i], num_t[k])
#         #         avg_time += t1[2]
#         #         avg_value += t1[1][1]
#         #     end
#         #     println("n= ", sizes[i], " p= ", p[j], " t= ", num_t[k])
#         #     println("avg_time: ", avg_time/10)
#         #     println("avg_cut:", avg_value/10)
#         # end
#     end
# end

"""
construct and test random graphs here
"""
for i =1:5
    for j = 1:4

        adj = []
        A = zeros(sizes[i],sizes[i])
        G = erdos_renyi(sizes[i], p[j])
        open(filename, "w") do file
            println(file, "$(sizes[i])")
            for k=1:sizes[i]
                nk = neighbors(G,k)
                push!(adj, nk)
                for l in nk
                    A[k,l] = 1
                    if l > k
                        println(file, "$k $l 1")
                    end
                end
            end
        end

       #  iters = [1,2,5,10,20,35,50]
       #  for k in iters
       #      t11 = @timed best_trevisan(A, adj, sizes[i], k)
       #      # push!(maxtrev_times, t11[2])
       #      # push!(maxtrev_cut_vals, t11[1][1])
       #      println("trevisan cut value with ", k, " iters: ", t11[1][1])
       #      println("trevisan time: ", t11[2])
       #      println()
       # end

        t0 = @timed simple_spectral(A,sizes[i],adj)
        # push!(spectral_times, t0[2])
        # push!(spectral_cut_vals, t0[1])
        println("simple spectral cut value: ", t0[1])
        println("simple spectral time: ", t0[2])

        t4 = @timed split_cuts(A,sizes[i],adj)
        # push!(split_times, t4[2])
        # push!(split_cut_vals, t4[1])
        println("split cuts cut value: ", t4[1])
        println("split cuts time: ", t4[2])


        t1 = @timed best_trevisan(A, adj, sizes[i], 5)
        # push!(trev_times, t1[2])
        # push!(trev_cut_vals, t1[1][1])
        println("trevisan cut value: ", t1[1][1])
        println("trevisan time: ", t1[2])

        t2 = @timed solve_max_cut_sdp(sizes[i], A, adj)
        # push!(sdp_times, t2[2])
        # push!(sdp_cut_vals, t2[1][1])
        println("SDP cut value: ", t2[1][1])
        println("SDP time: ", t2[2])

        t3 = @timed greedy_max_cut(A, collect(1:sizes[i]), adj, sizes[i])
        # push!(greedy_times, t3[2])
        # push!(greedy_cut_vals, t3[1][1])
        println("greedy cut value: ", t3[1][1])
        println("greedy time: ", t3[2])

    end
end


"""
test the TSPLIB graphs here
"""
# spectral_times = []
# spectral_cut_vals = []
# trev_times = []
# trev_cut_vals = []
# sdp_times = []
# sdp_cut_vals = []
# greedy_times = []
# greedy_cut_vals = []
#
for i=1:length(files)
    println(files[i])
    tsp = readTSPLIB(files[i])

    A = tsp.weights
    size = tsp.dimension
    adj = Array[]
    for s = 1:size
        push!(adj,filter(x -> x != s, collect(1:size)))
    end

    t0 = @timed simple_spectral(A,size,adj)
    # push!(spectral_times, t0[2])
    # push!(spectral_cut_vals, t0[1])
    println("simple spectral cut value: ", t0[1])
    println("simple spectral time: ", t0[2])

    t1 = @timed best_trevisan(A, adj, size, 5)
    # push!(trev_times, t1[2])
    # push!(trev_cut_vals, t1[1][1])
    println("trevisan cut value: ", t1[1][1])
    println("trevisan time: ", t1[2])

    t2 = @timed solve_max_cut_sdp(size, A, adj)
    # push!(sdp_times, t2[2])
    # push!(sdp_cut_vals, t2[1][1])
    println("SDP cut value: ", t2[1][1])
    println("SDP time: ", t2[2])

    t3 = @timed greedy_max_cut(A, collect(1:size), adj, size)
    # push!(greedy_times, t3[2])
    # push!(greedy_cut_vals, t3[1][1])
    println("greedy cut value: ", t3[1][1])
    println("greedy time: ", t3[2])

    t4 = @timed split_cuts(A, size ,adj)
    # push!(spectral_times, t0[2])
    # push!(spectral_cut_vals, t0[1])
    println("split cuts cut value: ", t4[1])
    println("split cuts time: ", t4[2])

    println()

end

"""
test the Network Repository graphs here
"""
# spectral_times = []
# spectral_cut_vals = []
# trev_times = []
# trev_cut_vals = []
# maxtrev_times = []
# maxtrev_cut_vals = []
# sdp_times = []
# sdp_cut_vals = []
# greedy_times = []
# greedy_cut_vals = []
file_names = []
file_sizes = []
open("files.txt") do f

  # read till end of file
  while ! eof(f)

     # read a new / next line for every iteration
     file = readline(f)

     adj = Array[]
     M = MatrixMarket.mmread("$file")
     size = M.m
     I,J,V = findnz(M)
     A = zeros(size,size)
     for s = 1:size
         push!(adj,Int[])
     end

     for k = 1:length(I)
         push!(adj[I[k]],J[k])
         A[I[k],J[k]] = V[k]
     end
     println(file)

     t0 = @timed simple_spectral(A,size,adj)
     # push!(spectral_times, t0[2])
     # push!(spectral_cut_vals, t0[1])
     println("simple spectral cut value: ", t0[1])
     println("simple spectral time: ", t0[2])

     t4 = @timed split_cuts(A,size,adj)
     # push!(spectral_times, t0[2])
     # push!(spectral_cut_vals, t0[1])
     println("split cuts cut value: ", t4[1])
     println("split cuts time: ", t4[2])

     t1 = @timed best_trevisan(A, adj, size, 5)
     # push!(trev_times, t1[2])
     # push!(trev_cut_vals, t1[1][1])
     println("trevisan cut value: ", t1[1][1])
     println("trevisan time: ", t1[2])

     t2 = @timed solve_max_cut_sdp(size, A, adj)
     push!(sdp_times, t2[2])
     push!(sdp_cut_vals, t2[1][1])
     println("SDP cut value: ", t2[1][1])
     println("SDP time: ", t2[2])

     t3 = @timed greedy_max_cut(A, collect(1:size), adj, size)
     # push!(greedy_times, t3[2])
     # push!(greedy_cut_vals, t3[1][1])
     println("greedy cut value: ", t3[1][1])
     println("greedy time: ", t3[2])

    println()
  end

end
