# fig6.jl
# Compare predictions against individual-based simulations for Human and
# Drosophila linkage maps.
# Arthur Zwaenepoel (2024)
using Pkg; Pkg.activate("/home/arthur_z/dev/Sewall")
using Distributed, Serialization
addprocs(6)
@everywhere using Sewall, Random
rng = Random.default_rng()
using Plots, PlotThemes, ColorSchemes; theme(:hokusai)

# Get the linkage maps and random architectures 
@everywhere begin
    Ls = 1.0          # total strength of selection
    L  = 100          # number of loci
    s  = Ls/L         # (haploid) selection coefficient
    u  = s/100        # per-locus mutation rate
    Ns = 10.0         # strength of selection relative to drif t
    Ne = Ns/s         # effective population size
    k  = 1            # number of diploids per haploid individual 
    N  = Ne2N(Ne, k)  # number of haploid individuals
    Ad = Sewall.randarch([Locus(-s, 0., 0., u, u) for _=1:L], flymap2)
    Ah = Sewall.randarch([Locus(-s, 0., 0., u, u) for _=1:L], humanmap2)
end

# Do the individual-based simulations
# Define the m/s values that we will examine for human and Drosophila
# respectively.
mss = [0.25:0.25:2.25, 0.25:0.25:3.0]   
res = map(zip([Ah, Ad], mss)) do (A, mss_)
    pmap(mss_) do ms
        seed  = rand(1:100)
        rng   = Random.seed!(seed)
        m     = ms*s
        D     = Deme(N=N, k=k, A=A)
        y     = ones(Bool, L)
        ml    = FixedMainland(y, (y, y))
        M     = MainlandIsland(D, m, 0.0, ml)
        P     = initpop(rng, M, zeros(Bool, L)) 
        P, xs = simulate!(rng, M, P, 120N, (P,_,_)->mean(P.haploids, dims=1))
        @info "Done $ms"
        seed, A, ms, vcat(xs[20N+2:end]...)
    end
end

# Get the analytical predictions
Qs = map(zip([Ah, Ad], mss)) do (A, mss_)
    m1, m2 = extrema(mss_)
    mssb = m1:0.01:m2
    qs = map(mssb) do ms
        @info ms
        m  = ms*s
        D  = Deme(N=N, k=k, A=A)
        y  = ones(Bool, L)
        ml = FixedMainland(y, (y, y))
        M  = MainlandIsland(D, m, 0.0, ml)
        Eq = mean.(Sewall.equilibrium(M))
    end
    Qs = hcat(qs...)
end

# Save the data
data = map(zip([Ah, Ad], mss, res, Qs)) do (A, mss_, sim, Q)
    m1, m2 = extrema(mss_)
    mssb = m1:0.01:m2
    ((mss=mss_, mssb=mssb, A=A, Ls=Ls, L=L, Ns=Ns, k=k, u=u), sim, Q)
end
serialize("data/linkage.jls", data)

# Load the results
(stuffh, resh, Qsh) = deserialize("data/linkage.jls")[1]
(stuffd, resd, Qsd) = deserialize("data/linkage.jls")[2]

# Make the plot
Pd = map(1:2) do i
    X     = [stuffd, stuffh][i]
    xmin  = [-0.02, -0.02][i]
    ytcks = [0.0:0.2:1, 0:0.2:1][i] 
    res   = [resd, resh][i]
    Qs    = [Qsd, Qsh][i] 
    slab  = ["\$Drosophila\$", "Human"][i]
    idx   = 15:20:X.L
    jdx   = [[2:2:length(res)-1 ; length(res)], 2:1:length(res)-3][i]
    s     = X.Ls/X.L
    rsm   = vec(s ./ mean(X.A.r[idx,:], dims=2))
    rsms  = map(x->@sprintf("%.2f", x), rsm)
    col   = collect(1:length(idx))'
    col2  = ColorSchemes.seaborn_colorblind[col]
    lab   = reshape(["locus$(@sprintf("%3d", i)) (\$s/\\bar{r}=$(rsms[j])\$)" for (j,i) in enumerate(idx)], 1, length(idx))
    # Comparison of all loci in the barrier
    P1 = plot(title="$slab, \$Ls=$(X.Ls), L=$(X.L), N_es=$(X.Ns)\$")
    map(res[jdx]) do (_,A,ms,x)
        i = findfirst(x->x==ms, X.mssb)
        zs = Qs[:,i]
        ys = 1 .- vec(mean(x, dims=1))
        scatter!(ys, zs, label="\$m/s=$(@sprintf("%.2f", ms))\$",
            legend=:topleft,
            xlabel="\$\\bar{p}\$ (simulation)", 
            ylabel="\$\\mathbb{E}[p]\$ (prediction)", ms=2)
    end
    plot!(x->x, size=(450,400), label="", color=:black, alpha=0.2, xlim=(xmin,1.001), ylim=(xmin,1.001), 
        yticks=ytcks, xticks=ytcks)
    # Couple of example loci across an `m/s` range
    P2   = plot(X.mssb, Qs[idx,:]', xlabel="\$m/s\$", color=col2, label=lab,
        legendfont=7, yticks=ytcks, ylim=(xmin,1))
    ys = map(xs->mean(xs[end][:,idx], dims=1), res)
    scatter!(X.mss, 1 .- vcat(ys...), 
        color=col2, label="", legend=:bottomleft,
        ylabel="\$\\mathbb{E}[p], \\bar{p}\$", ms=3)
    # The recombination map
    RR = [rij == 0.5 ? 0.5 : rij for rij in X.A.R] 
    SR = log2.(s ./ RR)
    mv = maximum(filter(!isnan, abs.(SR)))
    P3 = heatmap(SR, colormap=:vik, title="\$\\log_{2}s/r\$", clims=(-mv, mv),
        xlabel="locus", ylabel="locus", xlim=(-Inf,Inf), ylim=(-Inf,Inf))
    plot(P1, P2, P3, size=(800,220), layout=grid(1,3,widths=[0.28,0.28,0.44]), 
        bottom_margin=5Plots.mm, left_margin=3Plots.mm)
end
plot(reverse(Pd)..., layout=(2,1), size=(800,480))
