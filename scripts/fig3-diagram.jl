using Pkg
Pkg.activate(".")
using Graphs, Random, Arpack, Plots, LinearAlgebra, DataFrames, Statistics, CSV

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

αs = collect(0:0.1:4.1)
ks = collect(5:1.:999)

R = 100

df = DataFrame(α = Float64[], k_avg = Float64[], Δλ_avg = Float64[])

for α in αs
    for k in ks

        k_avg = 0.0
        Δλ_avg = 0.0

        for r in 1:R

            println("\u1b[1F")
            println(string(" α = ", α, " k = ", k, " R = ", r, "\u1b[1F"))

            # g = my_erdos_renyi(N, k)
            g = my_random_regular(N, k)
            
            try 

            fL = d_path_Laplacian(g, α)
            # fL = two_path_Laplacian(g, α)
            # fL = random_path_Laplacian(g, k, α)

            L = laplacian_matrix(g)

            fλ = eigs(fL, nev=3, which=:SM)[1][2]
            λ = eigs(L, nev=3, which=:SM)[1][2]

            Δλ = (fλ - λ)/N

            Δλ_avg += Δλ
            k_avg += mean(degree(g))

            catch
                continue
            end

        end

        Δλ_avg /= R
        k_avg /= R

        push!(df, [α, k_avg, Δλ_avg])
    end
end

CSV.write("fig3-diagram-RR.csv", df)