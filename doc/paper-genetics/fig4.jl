# fig4.jl
# Comparing single locus against multilocus predictions for heterogeneous
# genetic architectures.
# Arthur Zwaenepoel (2024)
using Sewall, Distributions, ThreadTools
using Plots, StatsBase, KernelDensity, StatsPlots, PlotThemes; theme(:hokusai)

# Function to obtain single locus predictions
function singlelocuseq(M, locus)
    A = Architecture([locus])
    D = Deme(N=M.deme.N, k=M.deme.k, A=A)
    M_ = MainlandIsland(D, M.mhap, M.mdip, M.mainland)
    mean(equilibrium(M_)[1])  # expected beneficial allele frequency on the island
end

# Function to generate multilocus and single locus predictions for random
# architectures drawn from a DFE model.
function randarches(dfe, L, N, k, u, m; n=100) 
    xs = [Float64[] for _=1:4]
    for i=1:n
        A  = Architecture([Locus(Sewall.randlocus(dfe)..., u, u) for i=1:L])
        ss = [l.s11 for l in A.loci]
        hs = [l.s01/l.s11 for l in A.loci]
        y  = ones(Bool, L)
        ml = FixedMainland(y, (y, y))
        M  = MainlandIsland(Deme(N=N, k=k, A=A), m, 0.0, ml)
        ds = Sewall.equilibrium(M)
        Ep = mean.(ds)
        ps = [singlelocuseq(M, locus) for locus in A.loci]
        push!(xs[1], ss...)
        push!(xs[2], hs...)
        push!(xs[3], Ep...)
        push!(xs[4], ps...)
    end
    return xs
end 

# Parameters
Lss = 0.5:0.5:1.5     # Ls values
mss = 0.05:0.15:0.5   # m/s values
Nes = 20              # strength of selection relative to drift 
s̄   = 0.01            # average selection coefficient
u   = 0.005*s̄         # mutation rate
Ne  = Nes/s̄           # effective population size 
k   = 1               # diploid individuals per haploid
N   = Ne2N(Ne, k)     # number of haploid individuals

# DFE model (independent selection and dominance coefficients in the diploid
# phase)
dfe = Sewall.IndependentDFE(Exponential(s̄), Beta())

# Obtain predictions for random architectures under the DFE
results = map(Lss) do Ls
    L = ceil(Int, Ls/s̄)
    nr = 150000 ÷ L
    @info Ls, L, nr
    tmap(mss) do ms
        xs = randarches(dfe, L, N, k, u, ms*s̄, n=nr)
    end
end

# Make the plot (panel A)
Ps = map(zip(results, Lss)) do (X, Ls)
    map(zip(X, mss)) do (xs, ms)
        idx = sample(1:length(xs[1]), 1000, replace=false)
        hs = xs[2][idx]
        zs = ColorSchemes.get.(Ref(ColorSchemes.viridis), hs)
        scatter(xs[4][idx], xs[3][idx], color=zs, alpha=0.7, ms=2,
            title="\$L\\bar{s} = $Ls, m/\\bar{s} = $ms\$",
            xlabel= Ls==1.5 ? "\$\\mathbb{E}[p_i]\$ single locus" : "",
            ylabel= ms==0.05 ? "\$\\mathbb{E}[p_i]\$ multilocus" : "",
        )
        plot!(x->x, color=:black, alpha=0.2, xlim=(0,1), ylim=(0,1))
    end
end |> x->vcat(x...) 
plot(Ps..., legend=false, size=(780,600))

# Function for fitting a KDE, assuming reflecting boundaries to correct for
# edge effects.
function reflectedkde(hs, ws, b1=0, b2=1)
    mx = maximum(hs)
    Hs = [b1 .- hs ; hs ; b2 .+ (mx .- hs)]
    Ws = [ws  ; ws ; ws]
    return kde(Hs, weights=Weights(Ws))
end

# Plot of the marginal densities at migration selection balance (B panel)
ps = map(zip(results, Lss)) do (X, Ls)
    xmx = 0.05; yh=(1-0.4,1+0.4)
    # plot the DFE of the genetic architecture (marginal distributions for s
    # and h)
    P1 = plot(0:0.001:xmx, s->pdf(dfe.sd, s), color=:black, xlim=(0,xmx),
        label="", fill=true, fillalpha=0.2)
    P2 = plot(0:0.01:1, h->Sewall.pdfh(dfe, h), color=:black, xlim=(0,1),
        ylim=yh, fill=true, fillalpha=0.2)
    # plot the DFE at migration selection balance (i.e. weighted by allele
    # frequency at equilibrium)
    map(enumerate(zip(X, mss))) do (i, (xs, m))
        ss  = abs.(xs[1])
        hs  = xs[2]
        ws  = xs[3]
        kd  = reflectedkde(hs, ws)
        D   = @sprintf "%.2f" mean(ws)
        ttl = "\$L\\bar{s}=$Ls\$"
        lab = @sprintf "%.2f" m
        p1  = density!(P1,
            ss, weights=Weights(ws), normalize=true, title=ttl,
            ylabel="\$f\$", label="\$m/\\bar{s} = $lab\$", xlabel="\$s_i\$",
            color=i, legend=Ls == 0.5 ? :topright : false)
        p2 = plot!(P2, 
            0:0.01:1, h->pdf(kd, h)*3, legend=false, ylabel="\$f\$",
            xlabel="\$h_i\$", color=i)
    end
    plot(P1, P2, layout=(2,1))
end |> x->hcat(x...)
ps = permutedims(ps)
plot(ps..., size=(200,850), layout=(3,1), left_margin=4Plots.mm)



