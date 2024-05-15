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
N  = 500
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
  N: Int64 500
  k: Int64 5
  Ne: Float64 454.5454545454545
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
 Wright{Float64}(Ne=454.5454545454545, A=0.03636363636363636, B=2.2321623989142325, sa=-0.008521220213390667, sb=-0.0008056822924104366, Z=80.80621883483394)
 Wright{Float64}(Ne=454.5454545454545, A=0.03636363636363636, B=2.215285608340794, sa=-0.0019059782417126887, sb=0.00046924455032741337, Z=27.172321361041828)
 Wright{Float64}(Ne=454.5454545454545, A=0.03636363636363636, B=2.2174481086583056, sa=-0.006693326373875632, sb=0.0027795370977858683, Z=31.439396859189273)
 Wright{Float64}(Ne=454.5454545454545, A=0.03636363636363636, B=2.2474791347605043, sa=-0.011030933973161334, sb=-0.0007019542143102373, Z=284.0852198940796)
 Wright{Float64}(Ne=454.5454545454545, A=0.03636363636363636, B=2.3018980811093273, sa=-0.022621727947580064, sb=-0.005251309607868718, Z=1.1107811376100693e7)
````

Expected equilibrium frequencies for the deleterious allelles:

````julia
mean.(ds[1:5])
````

````
5-element Vector{Float64}:
 0.4523421502616752
 0.026099277061638903
 0.0972126178097112
 0.6750689992008378
 0.8859013591444873
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
 0.18   0.996  0.176  0.182  0.126
 0.164  1.0    0.146  0.176  0.162
 0.144  1.0    0.15   0.184  0.142
 0.126  1.0    0.138  0.204  0.11
 0.146  1.0    0.098  0.216  0.124
 0.174  1.0    0.11   0.206  0.114
 0.176  0.998  0.09   0.218  0.12
 0.164  0.998  0.09   0.196  0.128
 0.182  0.998  0.092  0.204  0.12
 0.178  1.0    0.096  0.204  0.14
 0.198  1.0    0.086  0.214  0.114
````

Of course, one should do more/longer simulations.

If you have **any questions** on how to install, use, ... the package, *do
not hesitate to open an issue*, I'd be happy to help.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

