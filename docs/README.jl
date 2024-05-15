# # Sewall.jl

# This julia package implements:
#
# 1. Individual-based simulations of a multilocus, biallelic, mainland-island
#    or finite islands population genetics model for a
#    haplodiplontic/haploid/diploid life cycle.
# 2. Numerical tools to calculate effective migration rates and approximate
#    allele frequency distributions in the mainland-island/infinite islands
#    model using diffusion theory.
#
# To install the package:
#
# 1. Install `julia` from [https://julialang.org/](https://julialang.org/).
# 2. Open up a terminal and start a REPL by entering `julia` in the prompt.
# 3. Type `]` to enter the julia package manager.
# 4. Type `add https://github.com/arzwa/Sewall` to install the package.
#
# After executing these steps, one should be able to load the package using
# `using Sewall` in the julia REPL.

# Code for generating the figures in the main text of Zwaenepoel *et al.*
# (2024) can be found in `docs/figs-zwaenepoel2024`.

# ## Example

using Sewall, Distributions, Random; Random.seed!(12)

# This defines a genetic architecture with `L` biallelic unlinked loci with
# mean selection coefficient `s̄`, dominance coefficients distributed
# according to a Beta(2,2) distribution, and relative strength of haploid vs.
# diploid selection drawn from a uniform distribution. The `0->1` and `1->0`
# mutation rate at each locus is `u`.
function random_locus(s̄, u)
    s = rand(Exponential(s̄))
    h = rand(Beta(2,2))
    t = rand()
    return Locus(-s*(1-t), -s*h*t, -s*t, u, u)
end
L = 50
s̄ = 0.02
u = s̄*0.002
A = Architecture([random_locus(s̄, u) for i=1:L])
A[1:5]

# This defines a mainland-island model with `N` haploid individuals per
# generation, `N*k` diploid individuals per generation, genetic architecture `A`, haploid migration rate `m1`,
# diploid migration rate `m2` and mainland allele frequencies `y` 
N  = 550
k  = 5
D  = Deme(N=N, k=k, A=A)
m1 = s̄*0.4
m2 = 0.0
y  = ones(Bool, L)
ml = FixedMainland(y, (y, y))
M  = MainlandIsland(D, m1, m2, ml)

# Obtain the equilibrium allele frequency distributions (LE assumed)
ds = equilibrium(M)
ds[1:5]

# Expected equilibrium frequencies for the deleterious allelles:
mean.(ds[1:5])

# We can conduct an individual-based simulation for the same model
ngen = 1000
P = initpop(M, zeros(Bool, L))
P, xs = simulate!(M, P, ngen, (P,_,_)->mean(P.haploids, dims=1))
vcat(xs...)[end-10:end,1:5]

# Of course, one should do more/longer simulations.
#
# If you have **any questions** on how to install, use, ... the package, *do
# not hesitate to open an issue*, I'd be happy to help.

