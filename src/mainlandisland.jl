"""
    MainlandIsland
"""
@with_kw struct MainlandIsland{A,M,T}
    deme :: Deme{A}
    mhap :: T
    mdip :: T
    mainland :: M
    @assert nloci(mainland) == nloci(deme)
end

nloci(m::MainlandIsland) = length(m.deme)

"""
    FixedMainland{T}

A mainland which is fixed for a particular genotype.
"""
@with_kw struct FixedMainland{T}
    haploid :: Vector{T}
    diploid :: Tuple{Vector{T},Vector{T}}
    @assert length(haploid) == length(diploid[1]) == length(diploid[2])
end

diploidmigrant(_::AbstractRNG, m::FixedMainland) = m.diploid
haploidmigrant(_::AbstractRNG, m::FixedMainland) = m.haploid
nloci(m::FixedMainland) = length(m.haploid)

"""
    HWLEMainland

A mainland which is at HWLE with alelle frequencies `p`.
"""
struct HWLEMainland{T}
    p :: Vector{T}
end

nloci(m::HWLEMainland) = length(m.p)

function diploidmigrant(r::AbstractRNG, m::HWLEMainland)
    (haploidmigrant(r, m), haploidmigrant(r, m))
end

function haploidmigrant(r::AbstractRNG, m::HWLEMainland)
    rand(r, length(m.p)) .< m.p
end

function eqpdf(M::MainlandIsland)
    @unpack deme, mhap, mdip = M
    @unpack Ne, A = deme
    map(1:length(A)) do i
        @unpack u01, u10 = A[i]
        Wright(Ne, u01, mhap + mdip + u10, sasb(A[i])...)
    end
end

initpop(model, p0) = initpop(Random.default_rng(), model, p0)
initpop(rng::AbstractRNG, model::MainlandIsland, p0) = initpop(rng, model.deme, p0)

generation!(model, pop) = generation(Random.default_rng(), model, pop)
function generation!(rng::AbstractRNG, model::MainlandIsland, pop)
    @unpack deme, mhap, mdip, mainland = model
    haploidphase!(rng, model, pop)
    diploidphase!(rng, model, pop)
end


# XXX I guess it could be more efficient if we only track haploid genomes, keeping
# them twice in memory for copying (i.e. representing diploids as (i,j), where
# i, j are indices for the haploid genomes).
function haploidphase!(rng::AbstractRNG, model::MainlandIsland, pop)
    @unpack deme, mdip, mainland = model
    @unpack N, k = deme; Nk = N*k
    haploidfitness!(deme, pop)
    nm = rand(rng, Binomial(Nk, mdip)) 
    @views for i=1:nm
        mat, pat = diploidmigrant(rng, mainland)    
        pop.diploids[i   ,:] .= mat
        pop.diploids[Nk+i,:] .= pat
    end
    nk = Nk-nm
    idx = sample(rng, 1:N, Weights(pop.hfitness), 2nk, replace=true) 
    @views for j=1:nk
        pop.diploids[nm+j   ,:] .= pop.haploids[idx[j]   ,:]
        pop.diploids[Nk+nm+j,:] .= pop.haploids[idx[nk+j],:]
    end
end

function diploidphase!(rng::AbstractRNG, model::MainlandIsland, pop)
    @unpack deme, mhap, mainland = model
    @unpack N, k = deme; Nk = N*k
    diploidfitness!(deme, pop)
    nm = rand(rng, Binomial(N, mhap)) 
    @views for i=1:nm
        pop.haploids[i,:] .= haploidmigrant(rng, mainland)    
    end
    nk = N-nm
    idx = sample(rng, 1:Nk, Weights(pop.dfitness), nk, replace=true) 
    @views for j=1:nk
        mat = pop.diploids[idx[j]   ,:]
        pat = pop.diploids[Nk+idx[j],:]
        tgt = pop.haploids[nm+j     ,:]
        meiosis!(rng, tgt, deme.A, mat, pat) 
    end
    mutation!(rng, pop.haploids, deme.A)
end
