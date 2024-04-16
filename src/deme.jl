# Populations and demes
abstract type AbstractDeme end

"""
    Deme

Simple model for a haplodiplontic deme, with a discrete and synchronized
alternation of generations between N haploids and Nk diploids.
"""
@with_kw struct Deme{U<:Architecture} <: AbstractDeme
    N ::Int
    k ::Int
    Ne::Float64 = 1/(1/N + 1/(2N*k))
    A ::U  # genetic architecture
end

nloci(d::AbstractDeme) = length(d.A)
Ne2N(Ne, k) = ceil(Int, Ne/2k + Ne)

function eqpdf(d::Deme)
    @unpack Ne, A = d
    map(1:length(A)) do i
        @unpack u01, u10 = A[i]
        Wright(Ne, u01, u10, sasb(A[i])...)
    end
end

struct Population{T,V}
    haploids :: Matrix{T}
    diploids :: Matrix{T}
    hfitness :: Vector{V}
    dfitness :: Vector{V}
end

# We assume biallelic stuff, but we could if desired generalize by specializing
# function on the genetic architecture type (which should imply what represents
# an allele).
function initpop(rng::AbstractRNG, deme::AbstractDeme, p::AbstractVector)
    @unpack N, k = deme; L = length(p); Nk = N*k
    @assert L == nloci(deme) "Length of initial allele frequency vector does not match the genetic architecture"
    haploids = Matrix{Bool}(undef, N, L)
    diploids = Matrix{Bool}(undef, 2N*k, L)
    for i=1:N
        haploids[i,:] .= rand(rng, L) .< p
        for j=1:2k
            diploids[(i-1)*2k+j,:] .= rand(rng, L) .< p
        end
    end
    wh = Vector{Float64}(undef, N)
    wd = Vector{Float64}(undef, Nk)
    return Population(haploids, diploids, wh, wd)
end

function haploidfitness!(deme::Deme, pop::Population)
    pop.hfitness .= lognormalize(haploidfitness.(Ref(deme.A), eachrow(pop.haploids)))
end

function diploidfitness!(deme::Deme, pop::Population)
    @unpack N, k = deme; Nk = N*k
    @views for i=1:Nk
        pop.dfitness[i] = diploidfitness(deme.A, 
            pop.diploids[i,:], pop.diploids[Nk+i,:])
    end
    pop.dfitness .= lognormalize(pop.dfitness)
end

function haploidphase!(rng::AbstractRNG, deme::AbstractDeme, pop)
    @unpack N, k = deme; Nk = N*k
    haploidfitness!(deme, pop)
    idx = sample(rng, 1:N, Weights(pop.hfitness), 2Nk, replace=true)
    @views for i=1:2Nk
        pop.diploids[i,:] .= pop.haploids[idx[i],:]
    end
end

function diploidphase!(rng::AbstractRNG, deme::AbstractDeme, pop)
    @unpack N, k = deme; Nk = N*k
    diploidfitness!(deme, pop)
    idx = sample(rng, 1:Nk, Weights(pop.dfitness), N, replace=true)
    for i=1:N
        pat = @view pop.diploids[idx[i],:]
        mat = @view pop.diploids[idx[i]+Nk,:]
        meiosis!(rng, @view(pop.haploids[i,:]), deme.A, pat, mat)
    end
    mutation!(rng, pop.haploids, deme.A)
end

function generation!(rng::AbstractRNG, deme::AbstractDeme, pop)
    haploidphase!(rng, deme, pop)
    diploidphase!(rng, deme, pop)
end

# we can use the recombination probabilities directly, requires L rvs
# or we could use the map length, sample the number of crossovers from an
# appropriate Poisson distribution and then sample where they occur using a
# Multinomial?
@inline function meiosis!(rng::AbstractRNG, tgt, A::Architecture, x, y)
    pair = rand(rng) < 0.5 ? (x,y) : (y,x)
    k = 1
    for i=1:length(A)
        @inbounds tgt[i] = pair[k][i] 
        i == length(A) && break
        @inbounds k = rand(rng) < A.r[i] ? 1 + k%2 : k
    end
end

# This assumes pop contains haplotypes, i.e. length(pop) is the number of
# haplotypes. Note that each locus has its own forward and backward mutation
# rate.
function mutation!(rng::AbstractRNG, pop, arch)
    (N, L) = size(pop)
    ns = sum(pop, dims=1)
    for i in 1:L
        @unpack u01, u10 = arch[i]
        n10 = min(ns[i]  , rand(rng, Poisson(    ns[i]*u10)))
        n01 = min(N-ns[i], rand(rng, Poisson((N-ns[i])*u01)))
        idx = shuffle(1:N)
        k10 = k01 = 0
        k   = 1
        while true
            j = idx[k]
            if pop[j,i] && k10 < n10
                pop[j,i] = false
                k10 += 1
            elseif !pop[j,i] && k01 < n01
                pop[j,i] = true 
                k01 += 1
            end
            k += 1
            (k10 == n10 && k01 == n01) && break
        end
        # the above approach requires two Poisson rvs and a shuffle, and one
        # additional while loop which is somewhere in between O(nmut) and O(L).
        # I guess this is more efficient than sampling without replacement from
        # the two index sets? Should check it?
        #idx1 = findall(P.haploids[:,i])
        #idx0 = 1:L
        #idx = sample(rng, 1:N, nu, replace=false)  
        # sample without replacement? or with? shouldn't really matter
    end
end


