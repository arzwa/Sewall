# Locus definition
abstract type AbstractLocus end

"""
    Locus(s1, s01, s11)

Biallelic locus subject to selection in both the haploid and diploid phase of a
sexual life cycle.  Relative fitness of `0` allele is assumed to be 1 in both
haploids and diploid homozygotes.
"""
struct Locus{T} <: AbstractLocus
    s1  :: T  # haploid derived allele selection coefficient
    s01 :: T  # diploid heterozygous effect
    s11 :: T  # diploid homozygous effect
    u01 :: T  # mutation rate 0 -> 1
    u10 :: T  # mutation rate 1 -> 0
end

sasb(l::Locus) = (sa=l.s1 + l.s01, sb=l.s11 - 2l.s01)
sehe(l::Locus) = (se=2l.s1 + l.s11, 
    he=l.s01 + l.s11 == 0.0 ? 0.5 : (l.s1 + l.s01)/(2l.s1 + l.s11))

Base.show(io::IO, l::Locus) = write(io, "Locus$(string(l))")
Base.string(l::Locus) = @sprintf(
    "(%.3f, %.3f, %.3f, %.3f, %.3f)",
    l.s1, l.s01, l.s11, l.u01, l.u10)

# We are dealing with fitness on the log-scale. Note that we assume fitness is
# multiplicative across loci.
haploidfitness(l::Locus, g::Bool) = g ? l.s1 : 0.0
diploidfitness(l::Locus, h1::Bool, h2::Bool) = h1!=h2 ? l.s01 : (h1 ? l.s11 : 0.0)

# haplodiplontic mean fitness ∑ᵢpᵢeˢⁱ(∑ⱼpⱼeˢⁱʲ), reduces properly to diploid
# and haploid cases when the relevant selection coeffs are 0.
meanfitness(l::Locus, p, q=1-p) = p*(
    p+exp(l.s01)*q) + exp(l.s1)*q*(exp(l.s01)*p + exp(l.s11)*q)

# Genetic architecture
# TODO: should get rid of pairwise recombination rates -> they are not needed
# once we have the rate matrix...
"""
    Architecture

Genetic architecture, consists of a bunch of loci with certain effects, and a
linkage map (recombination rates). We assume `rrate[i]` gives the **probability
of a crossover** between locus `i` and locus `i+1`. 
"""
struct Architecture{T,V} <: AbstractVector{T}
    loci ::Vector{T}
    R    ::Matrix{V}  # recombination rates between all loci
end
function Architecture(loci::Vector{T}, r::Vector{V}) where {T<:AbstractLocus,V}
    @assert(length(r) == length(loci) - 1,
            "Recombination fraction vector does not match # of loci.")
    Architecture(loci, rrates(r))
end
Architecture(ls::Vector{Locus{T}}) where T = Architecture(ls, fill(0.5, length(ls)-1))
Architecture(l::Locus, L::Int) = Architecture([l for i=1:L], fill(0.5, L-1))
Architecture(l::Locus, L::Int, r) = Architecture([l for i=1:L], fill(r, L-1))

# AbstractVector methods
Base.size(A::Architecture) = (length(A.loci),)
Base.getindex(A::Architecture, i::UnitRange) = Architecture(A.loci[i], A.R[i,i]) 
Base.getindex(A::Architecture, i::Int) = A.loci[i]

haploidfitness(l::Architecture, g) = mapreduce(x->haploidfitness(x...), +, zip(l, g))
diploidfitness(l::Architecture, h1, h2) = mapreduce(x->diploidfitness(x...), +, zip(l, h1, h2))

#haploidfitness(A::Architecture, g) = haploidfitness(A.loci, g)
#diploidfitness(A::Architecture, g) = diploidfitness(A.loci, g)

#mappositions(A::Architecture) = [0.0; cumsum(invhaldane.(A.r))]
mappositions(A::Architecture) = [0.0 ; invhaldane.(A.R[2:end,1])]

#mappositions(r::Vector) = r == [0.5] ? [0.0] : [0.0; cumsum(invhaldane.(r))]

#function mappositions(A::Architecture)
#    idx = [1 ; findall(x->x==0.5, A.r)]
#    elems = [A.r[idx[i-1]:idx[i]-1] for i=2:length(idx)]
#    mappositions.(elems)
#end

function summarize(A::Architecture)
    L = length(A)
    w = countmap(A.loci)
    I = indexmap(A.loci)
    l = unique(A.loci)
    γ = [w[k] for k in l] 
    K = length(l)
    (loci=l, w=w, γ=γ, K=K, I=I, L=L) 
end

function lognormalize(l)
    x = exp.(l .- maximum(l))
    return x ./ sum(x)
end

function randarch(loci, M=humanmap)
    L = length(loci)
    _, _, rs = randloci(M, L)
    Architecture(loci, rs)
end

function randarch_neutral(loci, n, neutral_locus, M=GeneticMap(human_data()))
    all_loci, _, rs = randloci_neutral(M, loci, n, neutral_locus)
    Architecture(all_loci, rs)
end

