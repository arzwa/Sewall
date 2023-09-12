# Locus definition
"""
    Locus(s1, s01, s11)

Biallelic locus subject to selection in both the haploid and diploid phase of a
sexual life cycle.  Relative fitness of `0` allele is assumed to be 1 in both
haploids and diploid homozygotes.
"""
struct Locus{T}
    s1  :: T  # haploid derived allele selection coefficient
    s01 :: T  # diploid heterozygous effect
    s11 :: T  # diploid homozygous effect
    u01 :: T  # mutation rate 0 -> 1
    u10 :: T  # mutation rate 1 -> 0
end

sasb(l::Locus) = (sa=l.s1 + l.s01, sb=l.s11 - 2l.s01)
sehe(l::Locus) = (se=2l.s1 + l.s11, he=(l.s1 + l.s01)/(2l.s1 + l.s11))

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
"""
    Architecture

Genetic architecture, consists of a bunch of loci with certain effects, and a
linkage map (recombination rates). We assume `rrate[i]` gives the **probability
of a crossover** between locus `i` and locus `i+1`. 
"""
struct Architecture{T} <: AbstractVector{T}
    loci ::Vector{Locus{T}}
    r    ::Vector{T}  # recombination fractions between neighboring loci
    R    ::Matrix{T}  # recombination rates between all loci
    function Architecture(loci::Vector{Locus{T}}, r::Vector{T}) where T
        @assert(length(r) == length(loci) - 1,
                "Recombination fraction vector does not match # of loci.")
        new{T}(loci, r, rrates(r))
    end
end
Architecture(ls::Vector{Locus{T}}, r::Vector{V}) where {T,V} = Architecture(ls, map(T, r))
Architecture(ls::Vector{Locus}) = Architecture(ls, fill(0.5, length(ls)-1))
Architecture(l::Locus, L::Int) = Architecture([l for i=1:L], fill(0.5, L-1))
Architecture(l::Locus, L::Int, r) = Architecture([l for i=1:L], fill(r, L-1))

# AbstractVector methods
Base.size(A::Architecture) = (length(A.loci),)
Base.getindex(A::Architecture, i::UnitRange) = Architecture(A.loci[i], A.r[i[1:end-1]]) 
Base.getindex(A::Architecture, i::Int) = A.loci[i]
Base.vcat(A1::Architecture, A2::Architecture) = Architecture(vcat(A1.loci, A2.loci), vcat(A1.r, A2.r))

haploidfitness(l::Architecture, g) = mapreduce(x->haploidfitness(x...), +, zip(l, g))
diploidfitness(l::Architecture, h1, h2) = mapreduce(x->diploidfitness(x...), +, zip(l, h1, h2))


# Haldane's mapping function
haldane(y) = 0.5*(1-exp(-y))
invhaldane(x) = -log(1 - 2x)

# Compute recombination rates across a linear genome based on pairwise
# recombination fractions
rrates(xs) = hcat(rrates.(Ref(xs), 1:length(xs)+1)...)
function rrates(xs, j)
    ys = invhaldane.(xs)
    left = reverse(cumsum(ys[j-1:-1:1]))
    rght = cumsum(ys[j:end])
    [haldane.(left); NaN; haldane.(rght)]
end

#haploidfitness(A::Architecture, g) = haploidfitness(A.loci, g)
#diploidfitness(A::Architecture, g) = diploidfitness(A.loci, g)

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
