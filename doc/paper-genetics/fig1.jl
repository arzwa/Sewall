# fig1.jl 
# Evaluation of the approximation
# Arthur Zwaenepoel (2024)

# This is the code for figure 1, showing the fit of the multilocus diffusion
# approximation to data from indidvidual-based simulations.

# Load the required packages
using Pkg; Pkg.activate("/home/arthur_z/dev/Sewall")
using Distributed
addprocs(6)                       # add worker CPUs to run simulations in parallel
@everywhere using Sewall, Random  # load required packaged on all worker processes
using StatsBase
using Plots, PlotThemes           # packages for plotting

# Number of selected loci
L  = 50                           
# Overall strength of selection
Ls = 1.0
# Average strength of selection relative to drift $N_e\bar{s}$
Nes̄ = 10
# Average strength of selection at a locus
s̄  = Ls/L
# Number of diploid individuals per haploid individual
k  = 5
# Number of haploid indidviduals required to achieve the desired $N_e$
N  = Ne2N(Nes̄/s̄, k)                
# Get L random selection coefficients for the haploid phase 
s1s= rand(Exponential(s̄/2), L)
s1s .*= Ls/2sum(s1s)
# Get L random selection coefficients for the diploid phase 
s2s= rand(Exponential(s̄/2), L)
s2s .*= Ls/2sum(s2s)
# Get L random dominance coefficients for the diploid phase
hs = rand(L)
hs ./= 2mean(hs)
# Mutation rate
u  = s̄/100
# The loci that will constitute the genetic architecture
loci = [Locus(-s1, -s2*h, -s2, u, u) for (s1,s2,h) in zip(s1s,s2s,hs)]
# Values of m/s̄ for which we will run individual-based simulations
mss  = 0.05:0.1:2.5
# Number of generations we will simulate
ngen = 30N
# # ngen = 300N
# The generations we will sample to calculate expectations
keep = 5N:10:ngen
# # keep = 50N:10:ngen
@info length(keep)

# Do the simulations in parallel
res = pmap(mss) do ms
    seed = rand(1:100)
    rng = Random.seed!(seed)
    A  = Architecture(loci)
    m  = ms*s̄
    D  = Deme(N=N, k=k, A=A)
    y  = ones(Bool, L)
    ml = FixedMainland(y, (y, y))
    M  = MainlandIsland(D, m, 0.0, ml)
    P  = initpop(rng, M, zeros(Bool, L)) 
    ds = Sewall.equilibrium(M)
    P, xs = simulate!(rng, M, P, ngen, (P,_,_)->mean(P.haploids, dims=1))
    @info ms
    seed, A, ms, ds, xs 
end; 

# Calculate average allele frequencies
Xs = map(res) do (_, _, ms, ds, xs)
    XX = vcat(xs...)
    vec(mean(XX[keep,:], dims=1))
end

# Obtain the theoretical predictions
Ys = map(x->mean.(x[4]), res)

# Pick six random loci from the $m/\bar{s} = 0.45$ simulation to look at the
# allele frequency distribution.
k = findfirst(x->x>=0.45, mss)  # index of the simulation in the results
idx = rand(1:L, 6)

# Make the plot
P1 = plot(title="\$L=$L, \\bar{s}_1 = 0.01, \\bar{s}_{11} = 0.01, \\bar{h} = 0.5\$")
P2 = plot(xlabel="\$p\$", ylabel="\$\\log \\phi(p)\$", title="\$m/\\bar{s}=$(res[k][3])\$")
map(enumerate(idx)) do (j,i)
    lab = @sprintf "\$\\texttt{%.2f, %.2f, %.2f}\$" -100loci[i].s1 -100loci[i].s01 -100loci[i].s11
    plot!(P1, mss, getindex.(Ys, i), color=j, xlabel="\$m/\\bar{s}\$", label=lab,
        ylabel="\$\\mathbb{E}[p]\$", legend=false, xlim=extrema(mss))
    scatter!(P1, mss, 1 .- getindex.(Xs, i), color=j, label="", ms=2)
    d = res[k][4][i]
    X = getindex.(res[k][5], i)[keep]
    xs = 0.0001:0.003:0.9999
    plot!(P2, xs, map(p->logpdf(d, p), xs), color=j, legend=:bottomleft, label=lab, legendfontsize=6)
    hh = fit(Histogram, vec(X), 0:0.04:1,)
    hh = normalize(hh)
    es = hh.edges[1]
    ys = reverse(log.(hh.weights))
    xs = [es[i] + step(es) for i=1:length(ys)]
    scatter!(P2, xs, ys, color=j, ms=2, ylim=(-20,Inf), label="")
end
plot(P1, P2, size=(500,200), margin=3Plots.mm)



