# fig3.jl
# The effects of heterogeneity in selection coefficients (s)
# Arthur Zwaenepoel (2024)
using Pkg; Pkg.activate("/home/arthur_z/dev/Sewall")
using Distributed
addprocs(6)
@everywhere using Sewall, Random
rng = Random.default_rng()
using Plots, StatsPlots, PlotThemes; theme(:hokusai)

L    = 100
s̄    = 0.01
Nes  = 10.0
k    = 5
N    = Ne2N(Nes/s̄, k) 
u    = 0.005s̄
κs   = [8, 4, 2, 1, 1/2, 1/4]
mss  = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
nrep = 200

xs = pmap(mss) do ms
    @info ms
    ys = map(κs) do κ
        map(1:nrep) do _
            m  = ms*s̄
            dfe= Gamma(κ, s̄/κ)
            ss = rand(dfe, L)
            A  = Architecture([Locus(0., -ss[i]/2, -ss[i], u, u) for i=1:L])
            D  = Deme(N=N, k=k, A=A)
            y  = ones(Bool, L)
            ml = FixedMainland(y, (y, y))
            M  = MainlandIsland(D, m, 0.0, ml)
            ds = Sewall.equilibrium(M)
            Eq = mean(mean.(ds))
            gs = Sewall.eqgff(M, ds)
            Eq, gs
        end
    end
end;

# homogeneous
A  = Architecture([Locus(0., -s̄/2, -s̄, u, u) for _=1:L])
D  = Deme(N=N, k=k, A=A)
y  = ones(Bool, L)
ml = FixedMainland(y, (y, y))
ys = pmap(mss) do ms
    M  = MainlandIsland(D, ms*s̄, 0.0, ml)
    mean(Sewall.equilibrium(M)[1])
end

#al = ["8", "4", "2", "1", "\$\\frac{1}{2}\$", "\$\\frac{1}{4}\$"]
al = ["\$\\frac{1}{8}\$", "\$\\frac{1}{4}\$", "\$\\frac{1}{2}\$", "1", "2", "4"]
P1 = plot()
P2 = plot()
n1 = length(κs)
n2 = length(mss)
zs = map(enumerate(zip(mss,xs,ys))) do (j,(m,zs,y))
    ll = j == 1 ? "\$m/\\bar{s}=$m\$" : "\$$m\$"
    annotate!(P1, j*n1 + 0.6, 1.05, text(ll, 9, :right))
    plot!(P1, [0.5+6*(j-1),0.5+6*(j-1)+6], [y,y], color=:gray)
    zs = map(enumerate(zip(κs, zs))) do (i,(α,x_))
        qs = first.(x_)
        gs = map(mean, last.(x_)) ./ (m*s̄)
        # boxplot of mean barrier-wide 'multilocus predicted' differentiation across replicates
        boxplot!(P1, [6*(j-1)+i], qs, label="", fillalpha=0.9, lw=1, ms=2, whisker_width=0, color=i)
        boxplot!(P2, [6*(j-1)+i], gs, label="", fillalpha=0.9, lw=1, ms=2, whisker_width=0, color=i)
    end
end
vline!(P1, (n1+0.52):n1:n2*n1, color=:black, ylim=(0,1), label="")
vline!(P2, (n1+0.52):n1:n2*n1, color=:black, label="")
plot!(P1, 
    xticks=(1:n1*n2, repeat(fill("", 6), n2)), 
    xtickfont=7, 
    xlim=(0.5,n1*n2 + 0.5), 
    #xlabel="\$\\kappa^{-1} = \\mathrm{Var}[s]/\\bar{s}^2\$", 
    ylabel="\$\\bar{\\Delta}\$", size=(350,250),
    legend=false, top_margin=2Plots.mm, minorticks=false)
# This finishes plot 1
# Add the inset
xxs = 0.0001:0.0001:0.05
denss = mapreduce(κ->map(x->pdf(Gamma(κ, s̄/κ), x), xxs), hcat, κs)
labs  = reshape(["\$\\kappa = $κ\$" for κ in κs], (1,6))
plot!(xxs, denss, inset=bbox(0.04,0.1,0.24,0.34,:right), label=labs,
    legend=:topright, subplot=2, ylim=(0,200), tickfont=6, legendfont=6,
    ylabel="density", xlabel="\$s\$", guidefont=8, yticks=false)
# barrier strength
plot!(P2, xticks=(1:n1*n2, repeat(al, n2)), xtickfont=7, 
    xlim=(0.5,n1*n2 + 0.5), ylim=(0,1),
    xlabel="\$\\kappa^{-1} = \\mathrm{Var}[s]/\\bar{s}^2\$", 
    ylabel="\$\\bar{g} = \\overline{m_e}/m\$", size=(350,250),
    legend=false, top_margin=0Plots.mm, minorticks=false)
annotate!(P2, -2.8, 1.1, text("\$(\\mathrm{B})\$", 9))
annotate!(P1, -2.8, 1.05, text("\$(\\mathrm{A})\$", 9))
P12 = plot(P1, P2, layout=grid(2,1,heights=[0.7,0.3]), size=(450,350))
        
# differentiation across barrier
P4 = map(zip([0.1,0.4], ["C","D"])) do (ms,lab)
    P = plot(title="($lab) \$m/\\bar{s}=$ms\$", titlefont=9)
    Ps = map(κs) do κ
        m  = ms*s̄
        dfe= Gamma(κ, s̄/κ)
        ss = rand(dfe, L)
        A  = Architecture([Locus(0., -ss[i]/2, -ss[i], u, u) for i=1:L])
        D  = Deme(N=N, k=k, A=A)
        y  = ones(Bool, L)
        ml = FixedMainland(y, (y, y))
        M  = MainlandIsland(D, m, 0.0, ml)
        ds = Sewall.equilibrium(M)
        Eq = mean.(ds)
        plot!(P, sort(Eq, rev=true), line=:steppost, legend=:topright,
          ylabel="\$\\mathbb{E}[p]\$", xlabel="locus")
    end
    P
end |> x->plot(x..., legend=false, layout=(2,1))

plot(P12, P4, layout=grid(1,2,widths=[0.75,0.25]), size=(750,350),
    bottom_margin=2Plots.mm, left_margin=2Plots.mm)

