struct MainlandIsland{A,M,T}
    deme :: Deme{A}
    mhap :: T
    mdip :: T
    mainland :: M
end

struct FixedMainland{T}
    haploid :: Vector{T}
    diploid :: Tuple{Vector{T},Vector{T}}
end

sample(_::AbstractRNG, m::FixedMainland) = m.genotype

initpop(rng::AbstractRNG, model::MainlandIsland, p0) = initpop(rng, model.deme, p0)

function generation!(rng::AbstractRNG, model::MainlandIsland, pop)
    @unpack deme, mhap, mdip, mainland = model
    # haploid migration TODO
    haploidphase!(rng, deme, pop)
    # diploid migration TODO
    diploidphase!(rng, deme, pop)
end
