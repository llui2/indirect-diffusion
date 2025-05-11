using Pkg
Pkg.activate(".")
using Graphs, Random, Arpack, LinearAlgebra, DataFrames, Statistics, CSV

function d_path_Laplacian(g, α)
    n = nv(g)

    P_k = zeros(n,n)
    for i in 1:n
        P_k[i, :] = bellman_ford_shortest_paths(g, i).dists .^ (-α)
        P_k[i, i] = 0
    end

    Δ_k = Diagonal(P_k * ones(n))

    L_hat = Δ_k - P_k

    return L_hat
end

function two_path_Laplacian(g, α)
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

    Δ_k = zeros(n,n,2)
    for k in 1:2
        Δ_k[:,:,k] = Diagonal(transpose(transpose(ones(n)) * P_k[:,:,k]))
    end

    L_hat = zeros(n,n)
    for k in 1:2
        L_hat += float(k)^(-α) * ( Δ_k[:,:,k] - P_k[:,:,k] )
    end

    return L_hat
end

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

N = 1000
ks = collect((N-1)*log(N)/N:1.:999)

D = 10

df = DataFrame(k_avg = Float64[], Δλ_avg = Float64[], σ_Δλ = Float64[])

α = 0

for k in ks

    k_list = []
    Δλ_list = []

    for i in 1:D

        println("\u1b[1F")
        println(string(" k = ", k, " inst = ", i, "\u1b[1F"))

        G = my_erdos_renyi(N, k)
        # G = my_random_regular(N, k)

        L = two_path_Laplacian(G, α)
        # L = d_path_Laplacian(G, α)
        λ = eigs(L, nev=3, which=:SM)[1][2]

        nL = laplacian_matrix(G)
        nλ = eigs(nL, nev=3, which=:SM)[1][2]

        Δλ = (λ - nλ)/N

        push!(k_list, mean(degree(G)))       
        push!(Δλ_list, Δλ)

    end

    k_avg = mean(k_list)
    Δλ_avg = mean(Δλ_list)
    σ_Δλ = std(Δλ_list)

    push!(df, [k_avg, Δλ_avg, σ_Δλ])
end

CSV.write(string("fig-alphazero-ER.csv"), df)
