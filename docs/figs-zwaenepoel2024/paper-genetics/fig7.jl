# fig7.jl
# Simulate a genetic architecture and predict the Fₛₜ landscape, both
# accounting for LD and assuming LE among selected loci.
# Arthur Zwaenepoel (2024)
using Pkg; Pkg.activate("/home/arthur_z/dev/Sewall")
using Serialization, Parameters, Sewall, Random, StatsBase
rng = Random.default_rng()
using Plots, PlotThemes; theme(:hokusai)
fstr = "F_{\\mathrm{ST}}"

# Obtain the map locations for the loci.
# We will sample locations on the unit interval using a Dirichlet (this allows
# us to fine-tune to what extent the loci are more/less equally spaced), and
# then rescale to the entire map.
maplength = 10.0  # Map length in Morgans
α   = 10.0        # Dirichlet parameter (larger means more equally spaced)
xs  = cumsum(rand(Dirichlet(fill(α, T)))) .* maplength   # map locations
xs .-= xs[1]                       # start the map at zero
ds  = [xs[i] - xs[i-1] for i=2:T]  # pairwise distances  
rs  = Sewall.haldane.(ds)          # pairwise recombination rates

# Get the loci under selection.
T   = 1000  # total number of loci
L   = 100   # number of selected loci
Ls̄  = 0.7
s̄   = Ls̄/L
u   = s̄/100
m   = s̄*1.3
κ   = 1.0
dfe = Gamma(κ, s̄/κ)
selidx = sort(sample(1:T, L, replace=false))  # the indices for selected loci
neuidx = setdiff(1:T, selidx)                 # the indices for neutral loci
loci = [Locus(0.0, 0.0, 0.0, u, u) for i=1:T]
sss = rand(dfe, L)
sss .*= (s̄/mean(sss))
for (i,j) in enumerate(selidx)
    loci[j] = Locus(-sss[i], 0.0, 0.0, u, u)
end
A = Architecture(loci, rs)         # the genetic architecture

# Load the random architecture for Fig. 7, generated using the code above, from
# the CSV file.
using CSV, DataFrames
df = CSV.read("notebooks/paper-genetics/fig7-architecture.csv", DataFrame)
loci = [Locus(-s, 0.0, 0.0, u, u) for s in df[:,:s]]
A = Architecture(loci, df[1:end-1,:r])
xs = df[:,:x]
selidx = [i for i=1:T if A[i].s1 != 0.0] # the indices for selected loci
neuidx = setdiff(1:T, selidx)            # the indices for neutral loci

# Population size
Nes̄ = 25.0        
Ne  = Nes̄/s̄
k   = 1
N   = Ne2N(Ne, k)

# We assume mainland allele frequencies at approximate mutation-selection
# balance (neutral loci at 0.5).
y = [x.s1 == 0.0 ? 0.5 : max(min(0.9999, 1-(-u/x.s1)), 0.0001) for x in loci]
mainland = HWLEMainland(y)
D = Deme(N=N, k=k, A=A)
M = MainlandIsland(D, m, 0.0, mainland)

# Determine mₑ and expected Fₛₜ = (1 + 2Nₑmₑ)⁻¹
sld  = Sewall.singlelocuseq(M)
mld  = equilibrium(M) 
slme = Sewall.eqme(M, mean.(sld), Sewall.expectedpq.(sld))
mlme = Sewall.eqme(M)
mlFst = 1 ./ (2Ne .* mlme .+ 1)
slFst = 1 ./ (2Ne .* slme .+ 1) 

# Make the plot
ymax = maximum(mlFst)*1.5
xx = xs .* 100
tstr="(A) \$$fstr\$ landscape, \$L\\bar{s} = $Ls̄\$, \$m/\\bar{s} = $(round(M.mhap/s̄, digits=2)), N_e\\bar{s} = $(round(Ne*s̄)), \\kappa=$κ\$"
P1 = plot(ylabel="\$\\mathbb{E}[$fstr]\$", legend=:topleft,  ylim=(0,ymax), title=tstr)
plot!(xx[neuidx], mlFst[neuidx], label="accounting for LD")
plot!(xx[neuidx], slFst[neuidx], label="assuming LE")
ss = [-l.s1 for l in A] 
P0 = sticks(xx[ss .!= 0.0], ss[ss .!= 0.0], color=:black, legend=false, 
    ylabel="\$s\$", xlabel="map position (cM)", title="(B) selected sites")
P01 = plot(P1, P0, layout=grid(2,1,heights=[0.7,0.3]), size=(600,350),
    xlim=(0,1000), margin=2Plots.mm)

# Individual-based sims
# ngen = 10000
# reps = pmap(1:3) do _
#    MM = deepcopy(M)
#    init = rand(T) .< [x.s1 == 0.0 ? 0.5 : 0.0 for x in A]
#    P = initpop(rng, MM, init)     
#    ff(P, _, _) = vec(mean(P.diploids, dims=1))
#    pop, res = simulate!(rng, MM, P, ngen, ff)
#    X = hcat(res...)
#    (M, A, Ls, X)
# end
# zs = map(reps[1:3]) do rep
#    p = rep[end][neuidx,:]
#    fst = 1 .- p .* (1 .- p) ./ 0.25
# end
# Fst_ = mean(zs)
