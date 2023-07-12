include("trevisan_max_cut.jl")
include("official_SDP.jl")
include("helpers.jl")
include("split_cuts.jl")
include("burer.jl")
include("bmz_linesearch.jl")
using MatrixMarket
using LightGraphs
using Random
using SparseArrays
using TSPLIB
using PrettyTables

function burer_iter(it, A, size, adj, limit)
    cut = []
    cut_val = -Inf
    for i =1:it
        temp = burer_heuristic(A, size, adj, 10)
        if temp[1][1] > cut_val
            cut_val = temp[1]
            cut = temp[2]
        end
    end
    return cut_val, cut
end

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
    return local_search(A, adj_list, L, n, cut_value(A, L, adj_list))
end

"""
best_trevisan(A, adj, size, iter)
A : adjacency matrix of input graph
adj : adjacency list of input graph
size : number of vertices
iter : number of times to run trevisan's spectral algorithm

returns the best cut and value of 'iter' runs of trevisan's spectral algorithm
"""
function best_trevisan(A, adj, size, iter, limit)
    max = 0
    max_cut = []
    for i =1:iter
        cut, val = trevisan_max_cut(A, adj, collect(1:size), size, size)
        cut, val = local_search(A, adj, cut, size, val)
        if val > max
            max = val
            max_cut = cut
        end
    end
    return max, max_cut
end

#run these to prep the '@timed' function without having to process the function
@timed best_trevisan([0.0 1.0; 1.0 0.0], [[2],[1]], 2, 5, 1.0)
@timed solve_max_cut_sdp(2, [0.0 1.0; 1.0 0.0], [[2],[1]], false, 0, 10, true)
@timed simple_spectral([0.0 1.0; 1.0 0.0],2, [[2],[1]])
@timed greedy_max_cut([0.0 1.0; 1.0 0.0], collect(1:2), [[2],[1]], 2)
@timed split_cuts([0.0 1.0; 1.0 0.0], 2, [[2],[1]])
@timed burer_iter(1,[0.0 1.0; 1.0 0.0], 2, [[2],[1]], 10)

files = [:bayg29, :bays29, :berlin52, :bier127, :brazil58, :brg180, :ch130, :ch150, :d198, :eil101, :gr120, :gr137, :gr202, :gr96, :kroA100,:a280]

header = (["Graph\\Procedure", "Greedy", "Trevisan", "Simple Spectral", "Sweep Cuts", "SDP", "BMZ"])

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
burer_cut_vals = []
burer_times = []

file_names = []

"""
construct and test random graphs here
"""
#
open("random.txt") do f

  # read till end of file
  while ! eof(f)

     # read a new / next line for every iteration
    file = readline(f)
    println("$(file)")
    push!(file_names, "$(file)")
    open(file) do g
        graph = readline(g)
        sizes = split(graph, " ")
        n = parse(Int,sizes[1])
        m = parse(Int,sizes[2])
        A = zeros(n,n)
        adj = Array[]

        for i = 1:m
            temp = readline(g)
            edge = split(temp, " ")
            j = parse(Int,edge[1])
            k = parse(Int,edge[2])
            a = parse(Float64,edge[3])
            A[j,k] = 1
            A[k,j] = 1
        end

        for i = 1:n
            push!(adj, filter(j -> A[i,j] ==1, collect(1:n)))
        end
# #
        t0 = @timed simple_spectral(A,n,adj)
        push!(spectral_times, t0[2])
        push!(spectral_cut_vals, t0[1][2])
#
#
      avg_time = 0
      avg_cut = 0

      for p = 1:3
           t1 = @timed best_trevisan(A, adj, n, 5, 10)
           avg_time += t1[2]
           avg_cut += t1[1][1]
      end
      avg_time = avg_time/3
      avg_cut = avg_cut/3
      push!(trev_times, avg_time)
      push!(trev_cut_vals, avg_cut)
#
#       println("average trevisan time: ", avg_time)
#       println("average trevisan cut: ", avg_cut)
#
        avg_time = 0
        avg_cut = 0
#
        for p = 1:3
            t2 = @timed solve_max_cut_sdp(n, A, adj, false, 0, 10, true)
            avg_time += t2[2]
            avg_cut += t2[1][1]
            # println("sdp cut value: ", t2[1][1])
            # println("sdp time: ", t2[2])
        end
        avg_time = avg_time/3
        avg_cut = avg_cut/3

        push!(sdp_times, avg_time)
        push!(sdp_cut_vals, avg_cut)
#
#         println("average sdp time: ", avg_time)
#         println("average sdp cut: ", avg_cut)
#
        t3 = @timed greedy_max_cut(A, collect(1:n), adj, n)
        push!(greedy_times, t3[2])
        push!(greedy_cut_vals, t3[1][1])
        # println("greedy cut value: ", t3[1][1])
        # println("greedy time: ", t3[2])
#
        t4 = @timed split_cuts(A,n,adj)
        push!(split_times, t4[2])
        push!(split_cut_vals, t4[1])
#         println("split cuts cut value: ", t4[1])
#         println("split cuts time: ", t4[2])
#
#
# #
        avg_time = 0
        avg_cut = 0

        for p = 1:3

            t5 = @timed burer_iter(1, A, n, adj, 10)
            avg_time += t5[2]
            avg_cut += t5[1][1]

            # println("burer cut value: ", t5[1][1])
            # println("burer time: ", t5[2])
        end
          avg_time = avg_time/3
          avg_cut = avg_cut/3

          push!(burer_times, avg_time)
          push!(burer_cut_vals, avg_cut)
#         println("average burer time: ", avg_time/3)
#         println("average burer cut: ", avg_cut/3)
# #
      end
   end
end

random_time_data = hcat(file_names, greedy_times, trev_times, spectral_times, split_times, sdp_times, burer_times)
random_time_str = pretty_table(String, random_time_data; header = header)
output = open("output.txt", "w")
write(output, "Table 1")
write(output, random_time_str)


random_cut_data = hcat(file_names, greedy_cut_vals, trev_cut_vals, spectral_cut_vals, split_cut_vals, sdp_cut_vals, burer_cut_vals)
random_cut_str = pretty_table(String, random_cut_data; header = header)
write(output, "Table 2")
write(output, random_cut_str)

# println()





"""
test the TSPLIB graphs here
"""
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
burer_cut_vals = []
burer_times = []
file_names = []
#
for i=1:length(files)
    println(files[i])
    push!(file_names, files[i])
    tsp = readTSPLIB(files[i])

    A_mat = tsp.weights
    size = tsp.dimension
    adj = Array[]
    for s = 1:size
        push!(adj,filter(x -> x != s, collect(1:size)))
    end


    edgecount = (size*(size-1))/2
    open("$(files[i]).txt", "w") do file
        println(file, "$size $edgecount")
        for k=1:size
            for l = k+1:size
                println(file, "$k $l $(A_mat[k,l])")
            end
        end
    end
# #
# #
    t0 = @timed simple_spectral(A_mat,size,adj)
    push!(spectral_times, t0[2])
    push!(spectral_cut_vals, t0[1][2])
    # println("simple spectral cut value: ", t0[1][2])
    # println("simple spectral time: ", t0[2])

avg_time = 0
avg_cut = 0

for p = 1:3
    t1 = @timed best_trevisan(A_mat, adj, size, 5, 10)
    avg_time += t1[2]
    avg_cut += t1[1][1]

    # println("trevisan cut value: ", t1[1][1])
    # println("trevisan time: ", t1[2])
end
avg_time = avg_time/3
avg_cut = avg_cut/3
push!(trev_times, avg_time)
push!(trev_cut_vals, avg_cut)

# println("average trevisan time: ", avg_time)
# println("average trevisan cut: ", avg_cut)

avg_time = 0
avg_cut = 0

for p = 1:3
    t2 = @timed solve_max_cut_sdp(size, A_mat, adj, false, 0, 10, true)
    avg_time += t2[2]
    avg_cut += t2[1][1]

    # println("sdp cut value: ", t2[1][1])
    # println("sdp time: ", t2[2])
end
avg_time = avg_time/3
avg_cut = avg_cut/3
push!(sdp_times, avg_time)
push!(sdp_cut_vals, avg_cut)

# println("average sdp time: ", avg_time)
# println("average sdp cut: ", avg_cut)

# # #
    t3 = @timed greedy_max_cut(A_mat, collect(1:size), adj, size)
    push!(greedy_times, t3[2])
    push!(greedy_cut_vals, t3[1][1])
    # println("greedy cut value: ", t3[1][1])
    # println("greedy time: ", t3[2])

    t4 = @timed split_cuts(A_mat, size ,adj)
    push!(split_times, t4[2])
    push!(split_cut_vals, t4[1])
    # println("split cuts cut value: ", t4[1])
    # println("split cuts time: ", t4[2])
#


avg_time = 0
avg_cut = 0

for p = 1:3
    t5 = @timed burer_iter(1, A_mat, size, adj, 10)
    avg_time += t5[2]
    avg_cut += t5[1][1]

    # println("burer cut value: ", t5[1][1])
    # println("burer time: ", t5[2])
end
avg_time = avg_time/3
avg_cut = avg_cut/3
push!(burer_times, avg_time)
push!(burer_cut_vals, avg_cut)

# println("average burer time: ", avg_time)
# println("average burer cut: ", avg_cut)
#

end

tsp_time_data = hcat(file_names, greedy_times, trev_times, spectral_times, split_times, sdp_times, burer_times)
tsp_time_str = pretty_table(String, tsp_time_data; header = header)
write(output, "Table 3")
write(output, tsp_time_str)


tsp_cut_data = hcat(file_names, greedy_cut_vals, trev_cut_vals, spectral_cut_vals, split_cut_vals, sdp_cut_vals, burer_cut_vals)
tsp_cut_str = pretty_table(String, tsp_cut_data; header = header)
write(output, "Table 4")
write(output, tsp_cut_str)




"""
test the Network Repository graphs here
"""
spectral_times = []
spectral_cut_vals = []
trev_times = []
trev_cut_vals = []
sdp_times = []
sdp_cut_vals = []
greedy_times = []
greedy_cut_vals = []
burer_cut_vals = []
burer_times = []
split_times = []
split_cut_vals = []
file_names = []


open("files.txt") do f

  # read till end of file
  while ! eof(f)

     #read a new / next line for every iteration
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
         if I[k] != J[k] #ignore self-loops
             push!(adj[I[k]],J[k])
             if(A[I[k],J[k]] != 0)
                 println("duplicate")
             end
             A[I[k],J[k]] = V[k]
         end
     end
     println(file)
     push!(file_names, file)
#
     t0 = @timed simple_spectral(A,size,adj)
     push!(spectral_times, t0[2])
     push!(spectral_cut_vals, t0[1][2])
     # println("simple spectral cut value: ", t0[1][2])
     # println("simple spectral time: ", t0[2])


     avg_time = 0
     avg_cut = 0

     for p = 1:3
         t1 = @timed best_trevisan(A, adj, size, 5, 10)
         avg_time += t1[2]
         avg_cut += t1[1][1]

         # println("trevisan cut value: ", t1[1][1])
         # println("trevisan time: ", t1[2])
     end

     avg_time = avg_time/3
     avg_cut = avg_cut/3

     push!(trev_times, avg_time)
     push!(trev_cut_vals, avg_cut)

     # println("average trevisan time: ", avg_time)
     # println("average trevisan cut: ", avg_cut)



     avg_time = 0
     avg_cut = 0
     for p = 1:3
         t2 = @timed solve_max_cut_sdp(size, A, adj, false, 0, 10, true)
         avg_time += t2[2]
         avg_cut += t2[1][1]

         # println("sdp cut value: ", t2[1][1])
         # println("sdp time: ", t2[2])
     end
     avg_time = avg_time/3
     avg_cut = avg_cut/3

     push!(sdp_times, avg_time)
     push!(sdp_cut_vals, avg_cut)

     # println("average sdp time: ", avg_time)
     # println("average sdp cut: ", avg_cut)

      t3 = @timed greedy_max_cut(A, collect(1:size), adj, size)
      push!(greedy_times, t3[2])
      push!(greedy_cut_vals, t3[1][1])
      # println("greedy cut value: ", t3[1][1])
      # println("greedy time: ", t3[2])


     t4 = @timed split_cuts(A,size,adj)
     push!(split_times, t4[2])
     push!(split_cut_vals, t4[1])
     # println("split cuts cut value: ", t4[1])
     # println("split cuts time: ", t4[2])




     avg_time = 0
    avg_cut = 0
    for p = 1:3
        t5 = @timed burer_iter(1, A, size, adj, 10)
        avg_time += t5[2]
        avg_cut += t5[1][1]

        # println("burer cut value: ", t5[1][1])
        # println("burer time: ", t5[2])
    end
    avg_time = avg_time/3
    avg_cut = avg_cut/3

    push!(burer_times, avg_time)
    push!(burer_cut_vals, avg_cut)

    # println("average burer time: ", avg_time)
    # println("average burer cut: ", avg_cut)

  end

end

nr_time_data = hcat(file_names, greedy_times, trev_times, spectral_times, split_times, sdp_times, burer_times)
nr_time_str = pretty_table(String, nr_time_data; header = header)
write(output, "Table 5")
write(output, nr_time_str)


nr_cut_data = hcat(file_names, greedy_cut_vals, trev_cut_vals, spectral_cut_vals, split_cut_vals, sdp_cut_vals, burer_cut_vals)
nr_cut_str = pretty_table(String, nr_cut_data; header = header)
write(output, "Table 6")
write(output, nr_cut_str)

close(output)
