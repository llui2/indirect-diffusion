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

function my_scale_free(N::Int64, k::Float64)
    k = round(Int, k)
    m = Int(k * N)
    G = static_scale_free(N, m, 3)
    while !is_connected(G)
        G = static_scale_free(N, m, 3)
    end
    return G
end

function path_matrix(g)
    n = nv(g)
    
    A_k = zeros(n,n)
    for i in 1:n
        A_k[i, :] = bellman_ford_shortest_paths(g, i).dists
        A_k[i, i] = 0
    end

    d_max = 5
    P_k = zeros(n, n, d_max)
    
    # split each A_k
    for d in 1:d_max
        for i in 1:n
            for j in 1:n
                if A_k[i, j] == d
                    P_k[i, j, d] = 1
                end
            end
        end
    end

    return P_k
end

N = 1000
# ks = collect(2:1.:1000)
ks = exp10.(range(log10((N-1)*log(N)/N), log10(999), length=100))
R = 50

df = DataFrame(p1 = Float64[], p2 = Float64[], p3 = Float64[], p4 = Float64[], p5 = Float64[],
               sd_p1 = Float64[], sd_p2 = Float64[], sd_p3 = Float64[], sd_p4 = Float64[], sd_p5 = Float64[])

for k in ks

    p1s = []
    p2s = []
    p3s = []
    p4s = []
    p5s = []

    for r in 1:R
        
    println(string("k = ", k, ", r = ", r, "\u1b[1F"))

    g = my_erdos_renyi(N, k)
    # g = my_random_regular(N, k)
    # g = my_scale_free(N, k)

    P_k = path_matrix(g)

    A1 = P_k[:,:,1]
    A2 = P_k[:,:,2]
    A3 = P_k[:,:,3]
    A4 = P_k[:,:,4]
    A5 = P_k[:,:,5]

    g1 = Graph(A1)
    g2 = Graph(A2)
    g3 = Graph(A3)
    g4 = Graph(A4)
    g5 = Graph(A5)

    k1 = degree(g1)
    k2 = degree(g2)
    k3 = degree(g3)
    k4 = degree(g4)
    k5 = degree(g5)

    push!(p1s, mean(k1)/(N-1))
    push!(p2s, mean(k2)/(N-1))
    push!(p3s, mean(k3)/(N-1))
    push!(p4s, mean(k4)/(N-1))
    push!(p5s, mean(k5)/(N-1))

    end

    push!(df, [mean(p1s), mean(p2s), mean(p3s), mean(p4s), mean(p5s),
                std(p1s), std(p2s), std(p3s), std(p4s), std(p5s)])
end

CSV.write("fig2-probabilities-RR.csv", df)