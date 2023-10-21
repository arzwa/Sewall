"""
    MainlandIsland
"""
struct MainlandIsland{A,M,T}
    deme :: Deme{A}
    mhap :: T
    mdip :: T
    mainland :: M
end

nloci(m::MainlandIsland) = length(m.deme.A)

struct FixedMainland{T}
    haploid :: Vector{T}
    diploid :: Tuple{Vector{T},Vector{T}}
end

function eqpdf(M::MainlandIsland)
    @unpack deme, mhap, mdip = M
    @unpack Ne, A = deme
    map(1:length(A)) do i
        @unpack u01, u10 = A[i]
        Wright(Ne, u01, mhap + mdip + u10, sasb(A[i])...)
    end
end

diploidmigrant(_::AbstractRNG, m::FixedMainland) = m.diploid
haploidmigrant(_::AbstractRNG, m::FixedMainland) = m.haploid

initpop(model, p0) = initpop(Random.default_rng(), model, p0)
initpop(rng::AbstractRNG, model::MainlandIsland, p0) = initpop(rng, model.deme, p0)

generation!(model, pop) = generation(Random.default_rng(), model, pop)
function generation!(rng::AbstractRNG, model::MainlandIsland, pop)
    @unpack deme, mhap, mdip, mainland = model
    haploidphase!(rng, model, pop)
    diploidphase!(rng, model, pop)
end

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
