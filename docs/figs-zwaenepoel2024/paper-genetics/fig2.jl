# fig2.jl
# The effects of dominance
# Arthur Zwaenepoel (2024)

# Load the required packages
using Pkg; Pkg.activate("/home/arthur_z/dev/Sewall")
using Distributed
addprocs(6)                       # add worker CPUs to run simulations in parallel
@everywhere using Sewall, Random  # load required packaged on all worker processes
using StatsBase
using Plots, PlotThemes           # packages for plotting

# Selection coefficient
s   = 0.02 
# Overall strength of selection considered and number of loci
Ls  = [0.02, 0.5, 1.0, 1.5, 2.0] 
L   = ceil.(Int, Ls ./ s)
# Strength of selection relative to drift at a single locus
Nes = 20
Ne  = Nes / s
# Number of diploid individuals per haploid individual
k   = 2
# Number of haploid indidviduals required to achieve the desired $N_e$
N   = Ne2N(Nes/s, k)                
u   = s/100

# Values of m/s̄ for which we will run individual-based simulations
mss  = 0.05:0.1:1.2
# Number of generations we will simulate
NN   = ceil(Int, Ne)
ngen = 50NN
# # ngen = 300N
# The generations we will sample to calculate expectations
keep = 10NN:10:ngen
@info length(keep)

res = map([0.0, 0.5, 1.0]) do h
    map(L) do L_
        @info h, L_
        pmap(mss) do ms
            seed = rand(1:100)
            rng = Random.seed!(seed)
            A  = Architecture([Locus(0.0, -s*h, -s, u, u) for _=1:L_])
            D  = Deme(N=N, k=k, A=A)
            y  = ones(Bool, L_)
            ml = FixedMainland(y, (y, y))
            M  = MainlandIsland(D, ms*s, 0.0, ml)
            P  = initpop(rng, M, zeros(Bool, L_))
            P, xs = simulate!(rng, M, P, ngen, (P,_,_)->mean(P.haploids, dims=1))
            seed, M, ms, xs 
        end
    end
end

res2 = map([0.0, 0.5, 1.0]) do h
    map(L) do L_
        map(0:0.01:1.25) do ms
            A  = Architecture([Locus(0.0, -s*h, -s, u, u) for _=1:L_])
            D  = Deme(N=N, k=k, A=A)
            y  = ones(Bool, L_)
            ml = FixedMainland(y, (y, y))
            M  = MainlandIsland(D, ms*s, 0.0, ml)
            ms, mean(equilibrium(M)[1])
        end
    end
end

res3 = map(zip([0.0, 0.5, 1.0], res, res2)) do (h, Xs, Xs2)
    map(enumerate(zip(L, Xs, Xs2))) do (i, (L_, X, X2))
        ys = map(X) do (_, M, ms, xs)
            ms, 1-mean(vcat(xs...)[keep,:])
        end
        tx = getindex.(X2, 1)
        ty = getindex.(X2, 2)
        sx = getindex.(ys, 1)
        sy = getindex.(ys, 2)
        (L_*s, i, sx, sy, tx, ty)
    end
end

Ps = map(zip([0.0, 0.5, 1.0], res3)) do (h, Xs)
    P = plot(title="\$h=$h\$", legend=h == 0.5 ? :topright : false)
    for (lab, i, sx, sy, tx, ty) in Xs
        plot!(tx, ty, color=i, label="\$Ls = $lab\$")
        scatter!(sx, sy, color=i, label="", ms=3, 
            ylabel = h == 0.0 ? "\$\\mathbb{E}[p]\$" : "")
    end
    P
end
plot(Ps..., xlabel="\$m/s\$", layout=(1,3),
    size=(800,220), margin=6Plots.mm, ylim=(0,1),
    tickfont=10, titlefont=18, legendfont=10,
    guidefont=16, xlim=(0,1.25))


res_ = map([0.0, 0.5, 1.0]) do h
    L_ = L[end-1]
    map(0:0.01:1.25) do ms
        A  = Architecture([Locus(0.0, -s*h, -s, u, u) for _=1:L_])
        D  = Deme(N=2N, k=k, A=A)
        y  = ones(Bool, L_)
        ml = FixedMainland(y, (y, y))
        M  = MainlandIsland(D, ms*s, 0.0, ml)
        ms, mean(equilibrium(M)[1])
    end
end

