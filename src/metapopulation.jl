# Should do a new take on these: will probably be a lot more efficient to
# sample for instance haploid migrants directly from the diploid population of
# other localities of the previous generation proportional to local fitness
# instead of maintaining the conceptual distinction between selection
# and migration?

struct FiniteIslands{A,T}
    demes :: Vector{Deme{A}}
    mhap  :: Vector{T} # haploid migration rates
    mdip  :: Vector{T} # diploid migration rates
    Mhap  :: Matrix{T} # migration matrix haploids
    Mdip  :: Matrix{T} # migration matrix diploids
    function FiniteIslands(
            demes::Vector{Deme{A}},
            mhap ::Vector{T}, 
            mdip ::Vector{T}) where {T,A}
        d   = length(demes) 
        Ns  = [D.N for D in demes]
        Nks = [D.N*D.k for D in demes]
        Mhap = [mhap[i]/d .*  Ns[j] for i=1:d, j=1:d]
        Mdip = [mdip[i]/d .* Nks[j] for i=1:d, j=1:d]
        # there are two ways to see this: either m is the proportion of migrant
        # individuals, in which case the denominator should be d-1, or m is the
        # proportion of individuals drawn from the whole population (while 1-m
        # is the proportion of individuals drawn from the local population), in
        # which case d is the relevant denominator. The difference between the
        # two of course vanishes as d grows large... Li 1976 defines the
        # finite islands model in the latter way for instance. I think Himani
        # also formulates it like this, talking about a migrant pool.
        # In the latter case, m[i,i] = 1 - mᵢ(1 + 1/d)
        # XXX: actually it's a bit more tricky if we allow for unequal
        # popsizes, should double check...
        for i=1:d
            Mhap[i,i] += (1-mhap[i])*Ns[i]
            Mdip[i,i] += (1-mdip[i])*Nks[i]
            Mhap[i,:] = normalize(Mhap[i,:])
            Mdip[i,:] = normalize(Mdip[i,:])
        end
        new{A,T}(demes, mhap, mdip, Mhap, Mdip)
    end
end

normalize(xs) = xs ./ sum(xs)

function initpop(rng::AbstractRNG, model::FiniteIslands, p0)
    @unpack demes= model
    map(i->initpop(rng, demes[i], p0[i]), 1:length(demes))
end

function generation!(rng::AbstractRNG, model::FiniteIslands, pops)
    haploidphase!(rng, model, pops)
    diploidphase!(rng, model, pops)
end

# XXX Abstract out some patterns...
function haploidphase!(rng::AbstractRNG, model::FiniteIslands, pops)
    # confusingly, this involves the diploid migration phase...
    # i.e. it produces diploids, with selection and migration implemented in
    # one go
    @unpack demes, Mdip = model
    for (deme, pop) in zip(demes, pops)
        haploidfitness!(deme, pop)
    end
    for (d,(deme, pop)) in enumerate(zip(demes, pops))
        @unpack N, k = deme; Nk = N*k
        nmk = rand(rng, Multinomial(Nk, Mdip[d,:]))
        i = 1
        for (k,nk) in enumerate(nmk)  # nk from pop k
            # 2nk haplotypes that form the nk diploid migrants from pop k
            src = pops[k]
            NNk = demes[k].N
            idx = sample(rng, 1:NNk, Weights(src.hfitness), 2nk, replace=true)
            @views for j=1:nk
                pop.diploids[i   ,:] .= src.haploids[idx[j]   ,:]
                pop.diploids[Nk+i,:] .= src.haploids[idx[nk+j],:]
                i += 1
            end
        end
    end
end

function diploidphase!(rng::AbstractRNG, model::FiniteIslands, pops)
    @unpack demes, Mhap = model
    for (deme, pop) in zip(demes, pops)
        diploidfitness!(deme, pop)
    end
    for (d,(deme, pop)) in enumerate(zip(demes, pops))
        @unpack N, k = deme; Nk = N*k
        nmk = rand(rng, Multinomial(N, Mhap[d,:]))  # source demes for the migrants
        i = 1
        for (k,nk) in enumerate(nmk)  # nk from pop k
            # nk recombinant haplotypes that form the nk haploid migrants from pop k
            src = pops[k]
            Nkk = demes[k].N*demes[k].k
            idx = sample(rng, 1:Nkk, Weights(src.dfitness), nk, replace=true)
            @views for j=1:nk
                mat = src.diploids[idx[j]   ,:]
                pat = src.diploids[Nk+idx[j],:]
                tgt = pop.haploids[i,:]
                meiosis!(rng, tgt, demes[k].A, mat, pat) 
                i += 1
            end
        end
        mutation!(rng, pop.haploids, deme.A)
    end
end
