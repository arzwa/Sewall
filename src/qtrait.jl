"""
    QLocus
    
Quantitative trait locus for a biphasic life cycle, assuming the trait 'exists'
in both phases.
"""
struct QLocus{T} <: AbstractLocus
    a1::T  # additive effect in haploid phase
    a2::T  # homozygous additive effect in diploid phase
    d ::T  # dominance coefficient (∈ ℝ, i.e. `(1+d)a2/2` is the heterozygous effect) 
    u01::T
    u10::T
end

haploideffect(l::QLocus, g::Bool) = g ? l.a1 : 0.0
diploideffect(l::QLocus, h1::Bool, h2::Bool) = h1!=h2 ? l.a2*(1+l.d)/2 : (h1 ? l.a2 : 0.0)

haploidtrait(l::Architecture, g) = mapreduce(x->haploideffect(x...), +, zip(l, g))
diploidtrait(l::Architecture, h1, h2) = mapreduce(x->diploideffect(x...), +, zip(l, h1, h2))

"""
    QDeme

Simple model for a haplodiplontic deme, with a discrete and synchronized
alternation of generations between N haploids and Nk diploids, and a
quantitative trait under Gaussian stabilizing selection.
"""
@with_kw struct QDeme{U<:Architecture,T} <: AbstractDeme
    N ::Int
    k ::Int
    Ne::Float64 = 1/(1/N + 1/(2N*k))
    A ::U  # genetic architecture
    θh::T 
    θd::T
    Sh::T
    Sd::T 
end

haploidtraits(deme::QDeme, pop::Population) = haploidtrait.(
    Ref(deme.A), eachrow(pop.haploids))

function diploidtraits(deme::QDeme, pop::Population)
    @unpack N, k = deme; Nk = N*k
    map(i->diploidtrait(deme.A, pop.diploids[i,:], pop.diploids[Nk+i,:]), 1:Nk)
end

function haploidfitness!(deme::QDeme, pop::Population)
    @unpack θh, Sh = deme
    zh = haploidtraits(deme, pop)
    wh = -Sh .* ((zh .- θh) .^ 2)
    pop.hfitness .= lognormalize(wh)
end

function diploidfitness!(deme::QDeme, pop::Population)
    @unpack θd, Sd, N, k = deme; Nk = N*k
    @views for i=1:Nk
        zdi = diploidtrait(deme.A, pop.diploids[i,:], pop.diploids[Nk+i,:])
        pop.dfitness[i] = -Sd*(zdi - θd)^2
    end
    pop.dfitness .= lognormalize(pop.dfitness)
end


