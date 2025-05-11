using Pkg
Pkg.activate(".")
using Graphs, Random, Arpack, LinearAlgebra, DataFrames, Statistics, CSV

function my_erdos_renyi(N::Int64, k::Float64)
    p = k / (N - 1)
    G = erdos_renyi(N, p)
    while !is_connected(G)
        G = erdos_renyi(N, p)
    end
    return G
end

function my_random_regular(N::Int64, k::Float64)
    k = round(Int, k)
    G = random_regular_graph(N, k)
    while !is_connected(G)
        G = random_regular_graph(N, k)
    end
    return G
end

function path_matrix(g)
    n = nv(g)

    a1hop = adjacency_matrix(g)
    a2 = a1hop * a1hop
    a2 = a2 .> 0
    a2 = a2 .- a1hop
    a2 = max.(a2, 0)
    a2hop = a2 - Diagonal(Diagonal(a2))
    P_k = zeros(n,n,2)
    P_k[:,:,1] = a1hop
    P_k[:,:,2] = a2hop

    return P_k
end

N = 1000
k = 30.
D = 50

deg1_list = []
deg2_list = []

unit_vect = ones(N)

for i in 1:D

    println("\u1b[1F")
    println("i = ", i, "\u1b[1F")

    g = my_random_regular(N, k)
    # g = my_erdos_renyi(N, k)

    P_k = path_matrix(g)

    A1 = P_k[:,:,1]
    A2 = P_k[:,:,2]

    d1 = A1 * unit_vect
    d2 = A2 * unit_vect

    global deg1_list = vcat(deg1_list, d1)
    global deg2_list = vcat(deg2_list, d2)

end

df = DataFrame(deg1 = deg1_list, deg2 = deg2_list)

CSV.write("fig4-degreedist-RR.csv", df)