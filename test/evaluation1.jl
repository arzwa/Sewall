using Pkg; Pkg.activate("/home/arthur_z/dev/Sewall")
using Distributed
addprocs(6)
@everywhere using Sewall, Random
rng = Random.default_rng()
using Plots, PlotThemes; theme(:hokusai)

L  = 50
Ls = 1.0
Ns = 10
s̄  = Ls/L
k  = 5
N  = Ne2N(Ns/s̄, k)
s1s= rand(Exponential(s̄/2), L)
s2s= rand(Exponential(s̄/2), L)
s1s .*= Ls/2sum(s1s)
s2s .*= Ls/2sum(s2s)
hs = rand(L)
u  = s̄/100
loci = [Locus(-s1, -s2*h, -s2, u, u) for (s1,s2,h) in zip(s1s,s2s,hs)]

mss  = 0.05:0.1:2.5
ngen = 300N
keep = 50N:10:ngen
length(keep)

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
end

Xs = map(res) do (_, _, ms, ds, xs)
    XX = vcat(xs...)
    vec(mean(XX[keep,:], dims=1))
end

Ys = map(res) do (_, _, ms, ds, xs)
    mean.(ds)
end

k = 5
idx = rand(1:L, 6)
idx = [39,27,32,7,8,23]
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

