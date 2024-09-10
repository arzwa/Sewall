"""
    equilibrium(M::MainlandIsland, [p=ones(nloci(M))])

Determine the approximate marginal allele frequency distribution across the
L-locus architecture. Uses a fixed-point iteration initialized from `p`.  This
makes a lot of assumptions, among others loose linkage -- see Sachdeva (2022)
and Zwaenepoel et al. (2023).
"""
function equilibrium(M::MainlandIsland, p=ones(nloci(M)); tol=1e-9, kwargs...)
    migrant = GenePool(M.mainland)
    dists = pdfs(M, GenePool(p, p .* (1 .- p)), migrant)
    while true
        dists = pdfs(M, GenePool(p, expectedpq.(dists)), migrant)
        Ep = mean.(dists; kwargs...)
        norm(Ep .- p) < tol && return dists
        p = Ep
    end
end

"""
    eqgff(M::MainlandIsland [, ds::Vector{Wright}])

Return the equilibrium me values.
"""
eqgff(M::MainlandIsland) = eqgff(M, equilibrium(M))
eqgff(M::MainlandIsland, ds) = eqgff(M, mean.(ds), expectedpq.(ds))
function eqgff(M::MainlandIsland, Ep, Epq)
    migrant  = GenePool(M.mainland)
    resident = GenePool(Ep, Epq)
    @unpack mhap, mdip, deme = M
    @unpack A = deme
    gffs(A, mhap+mdip, resident, migrant)
end

function eqme(M::MainlandIsland, args...)
    @unpack mhap, mdip, deme = M
    gs = eqgff(M, args...)
    first.(gs) .* mhap .+ last.(gs) .* mdip
end

struct GenePool{T} # XXX rename?
    p ::Vector{T}
    pq::Vector{T}
end

function GenePool(m::FixedMainland)
    p = float.(m.haploid)  
    # XXX not sure what to do when haploid != diploid, but then this does not
    # make much sense for `FixedMainland`? or does it?
    pq = p .* (1 .- p)
    GenePool(p, pq)
end

function GenePool(m::HWLEMainland)
    p = m.p
    pq = p .* (1 .- p)
    GenePool(p, pq)
end

# XXX: this is not, in general, a good approximation.
function GenePool(m::Deme)
    ds = eqpdf(m)
    p  = 1 .- mean.(ds)
    pq = expectedpq.(ds)
    GenePool(p, pq)
end

function pdfs(M::MainlandIsland, focal, migrant)
    @unpack p, pq = migrant
    @unpack mhap, mdip, deme = M
    @unpack A, Ne = deme
    gs = gffs(A, mhap+mdip, focal, migrant)
    map(1:length(A)) do i
        s, h = sehe(A[i])
        me = gs[i][1]*mhap + gs[i][2]*mdip
        # XXX I guess when haploid and iploid mainland frequencies differ
        # we have to multiply mhap and mdip by the relevant frequencies
        # separately
        α = A[i].u10 + me*(1-p[i]) 
        β = A[i].u01 + me*p[i]
        Wright(Ne*s, Ne*α, Ne*β, h)
    end
end

gffs(A::Architecture, m, x, y) = map(j->gff(A, m, x, y, j), 1:length(A))

function gff(A::Architecture, m, focal::GenePool, migrant::GenePool, j::Int)
    @unpack p, pq = focal
    y, yz = migrant.p, migrant.pq
    x0 = x1 = 0.0
    for i=1:length(A)
        i == j && continue
        hap, dip = locuseffect(A[i], p[i], pq[i], y[i], yz[i], A.R[i,j], m)
        x0 += hap
        x1 += dip
    end
    g0 = exp(x0)
    g1 = g0*exp(x1)
    return g0, g1 
end

# The contribution of `locus` to the gff at a (possibly linked) locus with
# recombination probability `r`
function locuseffect(locus, p, pq, y, yz, r, m)
    @unpack s1, s01, s11 = locus
    sa, sb = sasb(locus)
    q = 1 - p
    z = 1 - y
    denom = m + r - sa*(p - q) - sb*(2pq - q)
    # should we do `denom = min(0.5,denom)`?
    hap = ( sa*(y - q) + sb*(pq - z*q)) / denom  # haploid migration
    dip = (s11*(y - q) + sb*(pq - yz ))  # diploid migration extra factor
    # note that the extra factor coming from diploid migrant selection does not
    # have anything in the denominator -> the denominator comes from the
    # series of backcrosses and is in `hap`!
    hap, dip
end

# (in)finite islands
function equilibrium(M::FiniteIslands; tol=1e-9, kwargs...)
    migrant = GenePool(M.mainland)
    dists = pdfs(M, GenePool(p, p .* (1 .- p)), migrant)
    while true
        dists = pdfs(M, GenePool(p, expectedpq.(dists)), migrant)
        Ep = mean.(dists)
        norm(Ep .- p) < tol && return dists
        p = Ep
    end
end

# Calculate an m_e profile by adding `n` neutral loci on the map and
# calculating local m_e.
# Need to deal with Inf map distances...
function meprofile(M::MainlandIsland, n; start=0.0, stop=0.0)
    # start with first locus at `start` M left from the current first locus
    # let last locus be `stop` M from last locus in architecture.
    @unpack A = M.deme
    RR, nxs, _, _, _ = get_profile_architecture(A, n, start, stop)
    ds = equilibrium(M)
    p  = mean.(ds)
    pq = expectedpq.(ds)
    # calculate
    migrant = GenePool(M.mainland)
    y, yz = migrant.p, migrant.pq
    m = M.mhap + M.mdip
    me = map(1:n) do j
        x0 = x1 = 0.0
        for i=1:length(A)
            hap, dip = locuseffect(A[i], p[i], pq[i], y[i], yz[i], RR[i,j], m)
            x0 += hap
            x1 += dip
        end
        g0 = exp(x0)
        g1 = g0*exp(x1)
        g0*M.mhap + g1*M.mdip 
    end
    nxs, me
end

"""
This calculates recombination rates for `n` hypothetical neutral loci scattered
on the genetic map defined by `A`. Returns both the L x n recombination rate
matrix (where the (i,j)th entry is the recombination rate for selected locus i
and neutral locus j) as well as the map positions of the neutral loci. This is
used for calculating m_e profiles.
"""
function get_profile_architecture(A::Architecture, n, start, stop)
    mappos = mappositions(A) 
    extra = range(start=0.0, stop=maximum(mappos) + start + stop, length=n)
    xs = [mappos .+ start; extra]
    o  = sortperm(xs)
    oi = invperm(o)
    nx = oi[length(mappos)+1:end]  # indices of the neutral mock loci 
    sx = oi[1:length(mappos)]  # indices of the selected loci
    xs = xs[o]
    # xs are the map locations of all loci
    R = ratematrix(xs)
    RR = [R[i,j] for i in sx, j in nx]
    return RR, extra, xs, sx, nx
end

function gff_aeschbacher(A::Architecture, n; start=0.0, stop=0.0)
    R, nxs, xs, sx, nx = get_profile_architecture(A, n, start, stop)
    gs = map(1:n) do j  # neutral locus j
        loggff = 0.0
        function _recurse1(a, i, j)
            (i > length(A) || sx[i] > nx[j]) && return 0.0
            b = _recurse1(a, i+1, j)
            loggff -= log(1 - A[i].s1/(R[i,j] + b))
            return a - b  # Aeschbacher expressions are for positive s...
        end
        function _recurse2(a, i, j)
            (i <= 0 || sx[i] < nx[j]) && return 0.0
            b = _recurse2(a, i-1, j)
            loggff -= log(1 - A[i].s1/(R[i,j] + b))
            return a - b
        end
        _recurse1(0.0, 1, j)
        _recurse2(0.0, length(A), j)
        loggff
    end
    nxs, exp.(gs)
end


