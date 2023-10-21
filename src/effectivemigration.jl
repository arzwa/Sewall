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
        Ep = mean.(dists)
        norm(Ep .- p) < tol && return dists
        p = Ep
    end
end

function eqgff(M::MainlandIsland)
    @unpack A = M.deme 
    ds = equilibrium(M)
    pm = float.(M.mainland.haploid)
    map(i->(ds[i].B - A[i].u01)/(2ds[i].Ne*pm[i]), 1:length(ds))
end

struct GenePool{T}
    p ::Vector{T}
    pq::Vector{T}
end

function GenePool(m::FixedMainland)
    p = float.(m.haploid)  # XXX not sure what to do when haploid != diploid 
    pq = p .* (1 .- p)
    GenePool(p, pq)
end

function pdfs(M::MainlandIsland, focal, migrant)
    @unpack p, pq = migrant
    @unpack mhap, mdip, deme = M
    @unpack A = deme
    gs = gffs(A, focal, migrant)
    map(1:length(A)) do i
        sa, sb = sasb(A[i])
        me = gs[i][1]*mhap + gs[i][2]*mdip
        # XXX I guess when haploid and iploid mainland frequencies differ
        # we have to multiply mhap and mdip by the relevant frequencies
        # separately
        α = A[i].u10 + me*(1-p[i]) 
        β = A[i].u01 + me*p[i]
        Wright(deme.Ne, α, β, sa, sb)
    end
end

gffs(A::Architecture, x, y) = map(j->gff(A, x, y, j), 1:length(A))

function gff(A::Architecture, focal::GenePool, migrant::GenePool, j::Int)
    @unpack p, pq = focal
    y, yz = migrant.p, migrant.pq
    x0 = x1 = 0.0
    for i=1:length(A)
        i == j && continue
        hap, dip = locuseffect(A[i], p[i], pq[i], y[i], yz[i], A.R[i,j])
        x0 += hap
        x1 += dip
    end
    g0 = exp(x0)
    g1 = g0*exp(x1)
    return g0, g1 
end

# The contribution of `locus` to the gff at a (possibly linked) locus with
# recombination probability `r`
function locuseffect(locus, p, pq, y, yz, r)
    @unpack s1, s01, s11 = locus
    sa, sb = sasb(locus)
    q = 1 - p; z = 1 - y
    hap = ( sa*(y - q) + sb*(pq - z*q)) / r  # haploid migration
    dip = (s11*(y - q) + sb*(pq - yz))  / r  # diploid migration extra factor
    hap, dip
end


