# Should do a new take on these: will probably be a lot more efficient to
# sample for instance haploid migrants directly from the diploid population of
# other localities of the previous generation proportional to local fitness
# instead of maintaining the conceptual distinction between selection
# and migration?

@with_kw struct FiniteIslands{A,T}
    demes :: Vector{Deme{A}}
    mhap  :: Vector{T}
    mdip  :: Vector{T}
    mhapp :: Vector{T} = normalize(mhap .* [d.N for d in demes])
    mdipp :: Vector{T} = normalize(mdip .* [d.N*d.k for d in demes])
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
    @unpack demes, mdip, mdipp = model
    for (deme, pop) in zip(demes, pops)
        haploidfitness!(deme, pop)
    end
    for (d,(deme, m, pop)) in enumerate(zip(demes, mdip, pops))
        @unpack N, k = deme; Nk = N*k
        nm = min(Nk, rand(rng, Poisson(Nk*m)))  # number of migrants
        nmk = rand(rng, Multinomial(nm, mdipp))  # source demes for the migrants
        nmk[d] += Nk - nm  # non-migrants
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
    @unpack demes, mhap, mhapp = model
    for (deme, pop) in zip(demes, pops)
        diploidfitness!(deme, pop)
    end
    for (d,(deme, m, pop)) in enumerate(zip(demes, mhap, pops))
        @unpack N, k = deme; Nk = N*k
        nm = min(N, rand(rng, Poisson(N*m)))  # number of migrants
        nmk = rand(rng, Multinomial(nm, mhapp))  # source demes for the migrants
        nmk[d] += N - nm
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
