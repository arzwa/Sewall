# fig5.jl
# Realized DFE at migration-selection balance for different DFE models for the
# underlying genetic architecture
# Arthur Zwaenepoel (2024)
using Sewall, Distributions, ThreadTools, QuadGK
using Plots, StatsBase, KernelDensity, StatsPlots, PlotThemes; theme(:hokusai)

mss = 0.05:0.15:0.8    # migration rates (relative to s̄)
nr  = 500              # number of replicate simulations
Nes = 20               # strength of selection relative to drift 
s̄   = 0.01             # average selection coefficient
u   = 0.005*s̄          # mutation rate
Ne  = Nes/s̄            # effective population size 
k   = 1                # diploid individuals per haploid
N   = Ne2N(Ne, k)      # number of haploid individuals
Ls  = 0.8              # total strength of selection
L   = ceil(Int, Ls/s̄)  # number of loci

# Define the different DFE models
κ = 1
λ = 1/s̄
dfe1 = Sewall.IndependentDFE(Gamma(κ, 1/λ), Beta(2,1))
dfe2 = Sewall.Logisticsbyh(Gamma(κ, 1/λ), (s̄/4, 0.5), (0.1135, 0.99), 1.)
dfe3 = Sewall.CKGamma(κ, λ, 1/3)

# Check that they have the same mean `h`
hm2 = quadgk(h->h*Sewall.pdfh(dfe2, h), 0, 1)[1]
hm3 = quadgk(h->h*Sewall.pdfh(dfe3, h), 0, 1)[1]
@info mean(dfe1.hd), hm2, hm3

# Function to generate multilocus predictions for random architectures drawn
# from a DFE model.
function randarches(dfe, L, N, k, u, m; n=100) 
    xs = [Float64[] for _=1:3]
    for i=1:n
        A  = Architecture([Locus(Sewall.randlocus(dfe)..., u, u) for i=1:L])
        ss = [l.s11 for l in A.loci]
        hs = [l.s01/l.s11 for l in A.loci]
        y  = ones(Bool, L)
        ml = FixedMainland(y, (y, y))
        M  = MainlandIsland(Deme(N=N, k=k, A=A), m, 0.0, ml)
        Ep = mean.(Sewall.equilibrium(M))
        push!(xs[1], ss...)
        push!(xs[2], hs...)
        push!(xs[3], Ep...)
    end
    return xs
end 

# Do the simulations
results = tmap(mss) do ms
    @info ms
    x1 = randarches(dfe1, L, N, k, u, ms*s̄, n=nr)
    x2 = randarches(dfe2, L, N, k, u, ms*s̄, n=nr)
    x3 = randarches(dfe3, L, N, k, u, ms*s̄, n=nr)
    (x1, x2, x3)
end

# Make the plots
Ps = map(enumerate(results)) do (j, x)
    pp = map(enumerate(zip([dfe1, dfe2, dfe3], x))) do (i, (dfe, y))
        ss, hs, ws = y
        D = sum(ws) / (nr*80)
        dd = @sprintf "%.2f" D
        ll = ["Independent", "Logistic", "CK94"][i]
        tt = "\$m/\\bar{s} = $(@sprintf "%.2f" mss[j])\$"
        xs = collect(zip(-ss, hs))
        ys = sample(xs, Weights(ws), 100000)
        kd = kde((first.(ys), last.(ys)), bandwidth=(0.005, 0.05))
        plot(kd, xlims=(0,0.06), ylims=(0,1), levels=7, colorbar=false,
             title  =   i == 1 ? tt : "",
             xticks = dfe == dfe3 ? (0:0.02:0.06) : false,
             xlabel = dfe == dfe3 ? "\$s\$" : "", 
             ylabel = j == 1 ? "$ll\n\$h\$" : "" )
        annotate!(0.058, 0.12, text("\$\\bar{\\Delta} = $dd\$", 8, :right))
    end
    plot(pp..., layout=(3,1), yticks=j == 1 ? (0:0.2:1) : false)
end 
plot(Ps..., layout=(1,length(mss)), 
    size=(750,350), left_margin=4Plots.mm,
    right_margin=-3Plots.mm, tickfont=7)

