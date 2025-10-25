module WrightDistribution

using QuadGK, Parameters, Reexport, LogExpFunctions
@reexport using Distributions
export Wright, sfs, expectedpq, Epq

"""
    Wright(Ne, α, β, sa, sb)

Wright's distribution for (haplo)diplontics, i.e.
```
f(p) ∝ p^(2Ne*α - 1) * (1-p)^(2Ne*β - 1) exp(-Ne(2sa*q + sb*q^2))
```
We assume `p` measures the `1` allele, so that `mean(d::Wright) = E[p]`.
"""
struct Wright{T} <: ContinuousUnivariateDistribution
    Ne :: T
    A  :: T 
    B  :: T
    sa :: T
    sb :: T
    Z  :: T
end

function Wright(Ne, α, β, sa, sb; kwargs...)
    A = 2Ne*α
    B = 2Ne*β
    Wright(promote(Ne, A, B, sa, sb, Zfun(Ne, A, B, sa, sb; kwargs...))...)
end

function Distributions.mean(d::Wright; kwargs...) 
	@unpack sa, sb, Ne, A, B, Z = d
    Yfun(Ne, A, B, sa, sb; kwargs...) / Z
end

function Distributions.var(d::Wright; kwargs...) 
    Ep = mean(d; kwargs...)
    Epq = expectedpq(d; kwargs...)
    Ep*(1-Ep) - Epq   
    # E[pq] = E[p] - E[p^2] <=> E[p^2] = E[p] - E[pq] 
    # => V[p] = E[p^2] - E[p]^2 
    #         = E[p] - E[pq] - E[p]^2 
    #         = E[p]E[q] - E[pq]  
end

Epq(d::Wright; kwargs...) = expectedpq(d; kwargs...)

function expectedpq(d::Wright; kwargs...)
	@unpack sa, sb, Ne, A, B, Z = d
    N, _ = quadgk(p -> p^A * (1-p)^B * Cfun(p, Ne, sa, sb), 0., 1.)
    return N/Z
end

function Distributions.pdf(d::Wright, p::Real) 
    p^(d.A-1) * (1-p)^(d.B-1) * Cfun(p, d.Ne, d.sa, d.sb) / d.Z
end

function Distributions.logpdf(d::Wright, p) 
    (d.A-1)*log(p) + (d.B-1)*log(1-p) - d.Ne*p*(2d.sa + d.sb*(2-p)) - log(d.Z)
end

function unnormpdf(d::Wright, p)
    p^(d.A-1) * (1-p)^(d.B-1) * Cfun(p, d.Ne, d.sa, d.sb)
end

function unnormlogpdf(d::Wright, p)
    (d.A-1)*log(p) + (d.B-1)*log(1-p) - d.Ne*p*(2d.sa + d.sb*(2-p))
end

function sfs(d::Wright, bins; kwargs...)
    ys = map(i->Zxfun(d, bins[i], bins[i+1]), 1:length(bins)-1)
    xs = bins[1:end-1] .+ step(bins)/2
    xs, ys
end

# Normalizing constant of Wright's distribution, using IBP
# could just use Z0fun/Z1fun/Zxfun below (but implemented this one first...)
function Zfun(N, A, B, sa, sb; kwargs...)
	c = 0.5^(A + B - 1) * exp(-N*(sa + 3sb/4)) * (1/A + 1/B)
	d, _ = quadgk(p->((    p^A)/A)*fpfun(p, B, N, sa, sb), 0, 0.5; kwargs...)
	e, _ = quadgk(p->(((1-p)^B)/B)*gpfun(p, A, N, sa, sb), 0.5, 1; kwargs...)
	c - d + e
end

function Yfun(N, A, B, sa, sb; kwargs...)
	C(p) = Cfun(p, N, sa, sb)
	L, _ = quadgk(p->p^A * (1-p)^(B-1) * C(p), 0., 0.5; kwargs...)
	R1 = (1/B)*(1/2)^(B+A)*C(0.5)
	R2, _ = quadgk(p->(((1-p)^B)/B)*gspfun(p, A, N, sa, sb), 0.5, 1; kwargs...)
	return L + R1 + R2
end

function Zxfun(d, x, y; kwargs...)
	@unpack sa, sb, Ne, A, B = d
    x ≈ 0 && y ≈ 1 && return d.Z
    x ≈ 0 && return Z0fun(d, y; kwargs...)
    y ≈ 1 && return Z1fun(d, x; kwargs...)
    I, _ = quadgk(p->pdf(d, p), x, y; kwargs...)
    return I
end

# 0 to x integral
function Z0fun(d, x; kwargs...)
	@unpack sa, sb, Ne, A, B = d
    a = (1/A)*(x)^A*(1-x)^(B-1)*Cfun(x, Ne, sa, sb)
	b, _ = quadgk(p->((p^A)/A)*fpfun(p, B, Ne, sa, sb), 0, x; kwargs...)
    return a - b
end

# x to 1 integral
function Z1fun(d, x; kwargs...)
	@unpack sa, sb, Ne, A, B = d
    #la = B * log(1-x) + (A-1)*log(x) + logCfun(x, Ne, sa, sb)
    #b, _ = quadgk(p->((1-p)^B)*gpfun(p, A, Ne, sa, sb), x, 1; kwargs...)
    #logsumexp(la, log(b))
    a = (1/B)*(1-x)^B*x^(A-1)*Cfun(x, Ne, sa, sb)
    b, _ = quadgk(p->(((1-p)^B)/B)*gpfun(p, A, Ne, sa, sb), x, 1; kwargs...)
    return a + b
end

Cfun(p, N, sa, sb) = exp(-N*p*(2sa + sb*(2-p)))
logCfun(p, N, sa, sb) = -N*p*(2sa + sb*(2-p))

function fpfun(p, B, N, sa, sb)
    -(1-p)^(B-1)*(N * (-p*sb + 2sa + sb*(2-p)) + (B-1)/(1-p)) * Cfun(p, N, sa, sb) 
end

function gpfun(p, A, N, sa, sb)
    p^(A-1)*(N*(p*sb - 2sa + sb*(p-2)) + (A-1)/p) * Cfun(p, N, sa, sb) 
end

function gspfun(p, A, N, sa, sb)
    p^A*(N*(p*sb - 2sa + sb*(p-2)) + A/p) * Cfun(p, N, sa, sb)
end


end # module WrightDistribution
