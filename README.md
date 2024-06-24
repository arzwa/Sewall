[![DOI](https://zenodo.org/badge/690597003.svg)](https://zenodo.org/doi/10.5281/zenodo.12515338)

# Sewall.jl

This julia package implements:

1. Individual-based simulations of a multilocus, biallelic, mainland-island
   or finite islands population genetics model for a
   haplodiplontic/haploid/diploid life cycle.
2. Numerical tools to calculate effective migration rates and approximate
   allele frequency distributions in the mainland-island/infinite islands
   model using diffusion theory.

To install the package:

1. Install `julia` from [https://julialang.org/](https://julialang.org/).
2. Open up a terminal and start a REPL by entering `julia` in the prompt.
3. Type `]` to enter the julia package manager.
4. Type `add https://github.com/arzwa/Sewall` to install the package.

After executing these steps, one should be able to load the package using
`using Sewall` in the julia REPL.

Code for generating the figures in the main text of Zwaenepoel *et al.*
(2024) can be found in `docs/figs-zwaenepoel2024`.

## Example

````julia
using Sewall, Distributions, Random; Random.seed!(12)
````

````
TaskLocalRNG()
````

This defines a genetic architecture with `L` biallelic unlinked loci with
mean selection coefficient `s̄`, dominance coefficients distributed
according to a Beta(2,2) distribution, and relative strength of haploid vs.
diploid selection drawn from a uniform distribution. The `0->1` and `1->0`
mutation rate at each locus is `u`.

````julia
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
````

````
5-element Architecture{Locus{Float64}, Float64}:
 Locus(-0.008, -0.001, -0.002, 0.000, 0.000)
 Locus(-0.001, -0.001, -0.002, 0.000, 0.000)
 Locus(-0.003, -0.004, -0.005, 0.000, 0.000)
 Locus(-0.011, -0.000, -0.001, 0.000, 0.000)
 Locus(-0.019, -0.004, -0.013, 0.000, 0.000)
````

This defines a mainland-island model with `N` haploid individuals per
generation, `N*k` diploid individuals per generation, genetic architecture `A`, haploid migration rate `m1`,
diploid migration rate `m2` and mainland allele frequencies `y`

````julia
N  = 550
k  = 5
D  = Deme(N=N, k=k, A=A)
m1 = s̄*0.4
m2 = 0.0
y  = ones(Bool, L)
ml = FixedMainland(y, (y, y))
M  = MainlandIsland(D, m1, m2, ml)
````

````
MainlandIsland{Architecture{Locus{Float64}, Float64}, FixedMainland{Bool}, Float64}(Deme{Architecture{Locus{Float64}, Float64}}
  N: Int64 550
  k: Int64 5
  Ne: Float64 500.0
  A: Architecture{Locus{Float64}, Float64}
, 0.008, 0.0, FixedMainland{Bool}(Bool[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], (Bool[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], Bool[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])))
````

Obtain the equilibrium allele frequency distributions (LE assumed)

````julia
ds = equilibrium(M)
ds[1:5]
````

````
5-element Vector{Wright{Float64}}:
 Wright{Float64}(Ne=500.0, A=0.04, B=2.4244157650758864, sa=-0.008521220213390667, sb=-0.0008056822924104366, Z=96.99890732528976)
 Wright{Float64}(Ne=500.0, A=0.04, B=2.403784755634715, sa=-0.0019059782417126887, sb=0.00046924455032741337, Z=24.59189118325925)
 Wright{Float64}(Ne=500.0, A=0.04, B=2.4064496849250308, sa=-0.006693326373875632, sb=0.0027795370977858683, Z=29.268602958523932)
 Wright{Float64}(Ne=500.0, A=0.04, B=2.440445276984162, sa=-0.011030933973161334, sb=-0.0007019542143102373, Z=435.0577321413137)
 Wright{Float64}(Ne=500.0, A=0.04, B=2.4979137568279763, sa=-0.022621727947580064, sb=-0.005251309607868718, Z=5.4348994895225875e7)
````

Expected equilibrium frequencies for the deleterious allelles:

````julia
mean.(ds[1:5])
````

````
5-element Vector{Float64}:
 0.5101708662142199
 0.027031413297720276
 0.10976065183051234
 0.7094968247768757
 0.8879364925871525
````

We can conduct an individual-based simulation for the same model

````julia
ngen = 1000
P = initpop(M, zeros(Bool, L))
P, xs = simulate!(M, P, ngen, (P,_,_)->mean(P.haploids, dims=1))
vcat(xs...)[end-10:end,1:5]
````

````
11×5 Matrix{Float64}:
 0.158182  0.963636  0.692727  0.152727   0.0727273
 0.149091  0.974545  0.689091  0.141818   0.0909091
 0.167273  0.981818  0.721818  0.165455   0.0927273
 0.205455  0.983636  0.716364  0.125455   0.0981818
 0.229091  0.985455  0.725455  0.123636   0.132727
 0.232727  0.983636  0.752727  0.0909091  0.103636
 0.247273  0.992727  0.727273  0.0927273  0.110909
 0.26      0.989091  0.738182  0.105455   0.129091
 0.256364  0.981818  0.792727  0.109091   0.14
 0.249091  0.983636  0.801818  0.12       0.138182
 0.254545  0.969091  0.785455  0.125455   0.145455
````

Of course, one should do more/longer simulations.

If you have **any questions** on how to install, use, ... the package, *do
not hesitate to open an issue*, I'd be happy to help.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

