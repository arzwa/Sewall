# stuff related to recombination maps
# Haldane's mapping function
haldane(d) = 0.5*(1-exp(-2d))     # distance -> recombination rate
invhaldane(r) = -0.5*log(1 - 2r)  # recombination rate -> distance

function maplength(A::Architecture)
    ml = 100sum(invhaldane.(A.r)) 
    @info "Map length = $ml cM (FYI human ~ 3800cM)"
end

# Compute recombination rates across a linear genome based on pairwise
# recombination fractions
rrates(xs) = hcat(rrates.(Ref(xs), 1:length(xs)+1)...)
function rrates(xs, j)
    ys = invhaldane.(xs)
    left = reverse(cumsum(ys[j-1:-1:1]))
    rght = cumsum(ys[j:end])
    [haldane.(left); NaN; haldane.(rght)]
end

# Drosophila data
drosophila_data_path = joinpath(dirname(pathof(Sewall)), 
    "..", "data/drosophila-rr-RRCv2.3/Comeron_tables/")

function drosophila_data()
    function _parse(x)
        x, d = split(x)
        parse(Int, x), parse(Float64, d)
    end
    rrs = map(readdir(drosophila_data_path, join=true)) do fname
        rr = map(_parse, readlines(fname))
        split(fname, "_")[end][1:end-4] => rr
    end |> Dict
    delete!(rrs, "chr4")  # there's no data for chr. 4 (tiny one)
end

# Human data
human_data_path = joinpath(dirname(pathof(Sewall)), 
    "..", "data/human-smoothed_map_b37_4/smoothed_map_b37_4/")

function human_data()
    function _parse(x)
        x, d = split(x, '\t')[[6,11]]
        parse(Int, x), parse(Float64, d)
    end
    rrs = map(readdir(human_data_path, join=true)) do fname
        rr = map(_parse, readlines(fname)[2:end])
        chr = parse(Int64, basename(fname)[4:end-6])
        chrname = @sprintf "chr%.2d" chr
        chrname => rr
    end |> Dict
    delete!(rrs, "chr23")  # consider autosomes only
    for (k,v) in rrs  # start the physical map where the genetic map starts (remove 'blind spots')
        minx = v[1][1]
        rrs[k] = [(x[1] - minx+1, x[2]) for x in v]
    end
    # the human maps record map units in cM vs. physical units (not map
    # units/physical unit)
    #                                     ↙ this is arbitrary 
    dd = Dict(k=>map(i->(v[i][1], max(1e-6, v[i+1][2] - v[i][2])), 
        1:length(v)-1) for (k,v) in rrs)
    return sort(dd)
end

# simpler genetic map (no local variation within chromosome)
human2_data_path = joinpath(dirname(pathof(Sewall)), 
    "..", "data/humanmap-matise/table1.csv")

function human2map()
    function _parse(x)
        x, d = split(x, ',')[[2,3]]
        parse(Float64, x), parse(Float64, d)
    end
    data = map(_parse, readlines(human2_data_path)[2:end])
    phys = round.(Int64, first.(data) .* 1e6)
    starts = cumsum([1; phys])
    starts = sort(vcat(starts, starts .-1))[2:end-1]
    dists = vcat([[last.(data)[i] ./ phys[i] ./ 100, Inf] for i=1:length(data)]...)
    chrs = [i:i+1 for i=1:2:2length(phys)] 
    GeneticMap(starts, dists, chrs)
end

function fly2map()
    @unpack starts, dists, chrs = flymap
    i = 1
    xs = map(chrs) do chr
        ds = [starts[i+1] - starts[i] for i=chr[1:end-1]]  # basepairs
        phys = sum(ds)  # chromosome length in Mb
        morgan = sum(dists[chr][1:end-1] .* ds) 
        ss = [starts[chr[1]], starts[chr[end]]]
        rd = [morgan/phys, Inf]
        cs = i:i+1 
        i += 2
        ss, rd, cs
    end
    GeneticMap(vcat(first.(xs)...), vcat(getindex.(xs, 2)...), last.(xs))
end


# We need some functions that will determine the rate of recombination between
# any two points on a given recombination map. This is fairly straightforward:
# calculate the distance in cM and use Haldane's mapping function to get
# recombination probabilities between all pairs of genes.
"""
    GeneticMap

Representation of a genetic map, with piecewise constant recombination rates.
`starts` record the starting points of a tract with a certain recombination
rate expressed in M/bp (i.e. `starts[i]` marks the starting point of a region
with recombination distance `dists[i]` M/bp).
Multiple chromosomes are represented by inserting an infinite distance between
two consecutive positions.
"""
struct GeneticMap{T,V,W}
    starts::Vector{T}  # bp
    dists ::Vector{V}  # M/bp
    chrs  ::Vector{W}  # index
end

Base.length(m::GeneticMap) = m.starts[end]

function GeneticMap(datadict)
    x = 0
    starts = Int64[]
    dists  = Float64[] 
    chrs   = UnitRange{Int64}[]
    x0 = 1
    for (_,v) in sort(datadict)
        # `k[1]` is the first base at which recombination rate `v` in cM/Mb
        # applies
        positions = x .+ first.(v)
        distances = last.(v) ./ 1e6 ./ 100  # XXX assumes cM/Mb
        x = maximum(positions) - 1
        # the last entry always has average recombination distance per bp 0,
        # here we insert Inf, to have recombination with 0.5 probability
        # between chromosomes.
        positions[end] -= 1
        distances[end] = Inf
        x1 = length(positions)
        starts = [starts ; positions]
        dists  = [dists  ; distances]
        push!(chrs, x0:x0+x1-1)
        x0 += x1
    end
    GeneticMap(starts, dists, chrs)
end

humanmap  = GeneticMap(human_data())
flymap    = GeneticMap(drosophila_data())
humanmap2 = human2map()
flymap2   = fly2map()

function distance(m::GeneticMap, x1, x2)
    @unpack starts, dists = m
    x1, x2 = min(x1,x2), max(x1,x2)
    physd = x2 - x1
    i1 = findfirst(i->starts[i] > x1, 1:length(starts))-1
    i2 = findlast( i->starts[i] < x2, 1:length(starts))
    # exception to the below
    i1 == i2 && return dists[i1]*physd
    #    physical distance to the end of the first part
    d1 = (starts[i1+1]-1) - x1
    #    physical distance after the beginning of the last part
    d2 = x2 - starts[i2]
    d  = d1*dists[i1] + d2*dists[i2]
    for j=i1+1:i2-1
        d += dists[j]*(starts[j+1] - starts[j-1])
    end
    return d
end

"""
    randloci(GeneticMap, L)

Pick L random loci along a given genetic map and get their pairwise
recombination rates.
"""
function randloci(m::GeneticMap, L) 
    locs = sort(sample(1:length(m), L, replace=false))
    ds = map(1:L-1) do i
        distance(m, locs[i], locs[i+1])
    end
    locs, ds, haldane.(ds)
end

#function equallyspacedloci(m::GeneticMap, L)
#end

"""
    randloci_neutral(GeneticMap, L, n)

Pick L random selected loci and n regularly spaced neutral loci along a given
genetic map and get their pairwise recombination rates.
"""
function randloci_neutral(m::GeneticMap, loci, n, neutral) 
    ml = length(m)
    L  = length(loci)
    step = ceil(Int, ml / n)
    neutral_loci = 1:step:ml 
    selected_loci = Int[]
    while length(selected_loci) < L 
        k = sample(1:ml)
        (k ∉ selected_loci && k ∉ neutral_loci) && push!(selected_loci, k)
    end
    all_idx = [selected_loci; neutral_loci]
    ds = map(1:L+n-1) do i
        distance(m, all_idx[i], all_idx[i+1])
    end
    o = sortperm(all_idx)
    all_loci = [loci ; [neutral for _=1:n]][o]
    all_loci, ds, haldane.(ds)
end

function maplength(g::GeneticMap)
    @unpack starts, dists, chrs = g
    map(chrs) do chr
        ds = [starts[i+1] - starts[i] for i=chr[1:end-1]]  # basepairs
        phys = sum(ds)  # chromosome length in Mb
        morgan = sum(dists[chr][1:end-1] .* ds) 
        (phys, phys/1e6, morgan*100)
    end
end

