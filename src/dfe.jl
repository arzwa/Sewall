# Distributions of fitness effects
struct IndependentDFE{T,V}
    sd::T
    hd::V
end

function randlocus(d::IndependentDFE)
    s = rand(d.sd)
    h = rand(d.hd)
    0.0, -s*h, -s
end

Distributions.pdf(d::IndependentDFE, l)    = pdf(d, -l.s11, l.s01/l.s11)
Distributions.pdf(d::IndependentDFE, s, h) = pdf(d.sd, s) * pdf(d.hd, h)
Distributions.logpdf(d::IndependentDFE, l)    = logpdf(d, -l.s11, l.s01/l.s11)
Distributions.logpdf(d::IndependentDFE, s, h) = logpdf(d.sd, s) + logpdf(d.hd, h)

pdfh(dfe::IndependentDFE, h) = pdf(dfe.hd, h)

# Logistic relationship
struct Logisticsbyh{T,V}
    sd::V
    a::T
    b::T
    σ::T
end

function Logisticsbyh(sd, sh1::T, sh2::T, σ) where T<:Tuple 
    a, b = _getab(sh1..., sh2...)
    Logisticsbyh(sd, a, b, σ)
end

function randlocus(dfe::Logisticsbyh)
    @unpack sd, a, b, σ = dfe
    s = rand(sd)
    h = _hfroms(s, a, b, σ)
    0.0, -s*h, -s
end

function _getab(s1, h1, s2, h2)
    b = (log(h2/(1-h2)) - log(h1/(1-h1)))/(log(s2) - log(s1))
    a = log(h1/(1-h1)) - b*log(s1)
    a, b
end

function _hfroms(s, a, b, σ)
    lh = a + b*log(s) + rand(Normal(0,σ))
    h  = 1/(1 + exp(-lh))
end

Distributions.pdf(d::Logisticsbyh, l) = pdf(d, -l.s11, l.s01/l.s11)
Distributions.pdf(d::Logisticsbyh, s, h) = 
    pdf(d.sd, s) * pdf(Normal(d.a + d.b*log(s), d.σ), logit(h)) / (h*(1-h)) 
Distributions.logpdf(d::Logisticsbyh, l) = logpdf(d, -l.s11, l.s01/l.s11)
Distributions.logpdf(d::Logisticsbyh, s, h) = 
    logpdf(d.sd, s) + logpdf(Normal(d.a + d.b*log(s), d.σ), logit(h)) - log(h) - log(1-h)

# p(h) = 1/(h(1-h))*p(lh) = 1/(h(1-h))∫p(lh|s)p(s)ds = 1/(h(1-h))∫N(lh|a + b*log(s),σ)p(s)ds
function pdfh(dfe::Logisticsbyh, h)
    @unpack sd, a, b, σ = dfe
    l = log(h/(1-h))
    I, _ = quadgk(s->pdf(sd, s)*pdf(Normal(a + b*log(s), σ), l), 0, 1.)
    return I/(h*(1-h))
end

# Caballero & Keightley 1994 (apparently, read in Zhang & Hill)
# Note that in Agrawal & Whitlock 2011, this model is criticized: "On the basis
# of a comparatively small number of transposable element insertions in
# Drosophila, Caballero and Keightley (1994) made some inferences about the
# joint distribution of h and s. They proposed that h, for a given s, was
# uniformly distributed between 0 and exp(-Ks). Our results differ from this
# model in two ways with respect to the variance in h.  First, our results in-
# dicate that, for a given s, values of h are not uniformly distributed but
# rather strongly skewed. Second, the Caballero and Keightley model predicts a
# decline in the variance in h with increasing s; we found no support for such
# a decline in the yeast data set."

# CKGamma
struct CKGamma{T,V<:Distribution}
    sd::V
    h̄::T
    K::T
end
    
CKGamma(κ, λ, h̄) = CKGamma(Gamma(κ,1/λ), h̄, λ*((2h̄)^(-1/κ) - 1))

function randlocus(d::CKGamma) 
    s = rand(d.sd)
    h = 1-rand(Uniform(0, exp(-d.K*s)))  
    # note we use 1-h, not h, to get a large-effects tend to be recessive
    # pattern (i.e. for large effect alleles, the wild type coming from the
    # mainland tends to be dominant, these large effect alleles tended to be
    # recessive on the mainland...)
    0.0, -s*h, -s
end

Distributions.pdf(d::CKGamma, l) = pdf(d, -l.s11, l.s01/l.s11)
Distributions.pdf(d::CKGamma, s, h) = 1-h > exp(-s*d.K) ? 0. : pdf(d.sd, s) * exp(s*d.K) 
Distributions.logpdf(d::CKGamma, l) = log(pdf(d, l))
Distributions.logpdf(d::CKGamma, s, h) = log(pdf(d, s, h))
 
function pdfh(d::CKGamma, h)
    β = 1/d.sd.θ
    α = d.sd.α
    K = d.K
    (β/(β - K))^α * cdf(Gamma(α, 1/(β - K)), -log(1-h)/K) 
end

#struct CKGamma2{T,V<:Distribution}
#    sd::V
#    h̄::T
#    K::T
#end
#    
#CKGamma2(κ, λ, h̄) = CKGamma2(Gamma(κ,1/λ), h̄, λ*((2h̄)^(-1/κ) - 1))
#
#function randlocus(d::CKGamma2) 
#    s = rand(d.sd)
#    h = rand(Uniform(0, exp(-d.K*s)))  
#    DipLocus(-s*h, -s)
#end
#
#Distributions.pdf(d::CKGamma2, l) = pdf(d, -l.s11, l.s01/l.s11)
#Distributions.pdf(d::CKGamma2, s, h) = h > exp(-s*d.K) ? 0. : pdf(d.sd, s) * exp(s*d.K) 
#Distributions.logpdf(d::CKGamma2, l) = log(pdf(d, l))
#Distributions.logpdf(d::CKGamma2, s, h) = log(pdf(d, s, h))
# 
#function pdfh(d::CKGamma2, h)
#    β = 1/d.sd.θ
#    α = d.sd.α
#    K = d.K
#    (β/(β - K))^α * cdf(Gamma(α, 1/(β - K)), -log(h)/K) 
#end
#
#struct AWGamma{T,V<:Distribution,W<:Distribution}
#    sd::V
#    hd::W
#    β1::T
#    β2::T
#end
#
#AWGamma(sd, δ, Vh, β1, β2) = AWGamma(sd, Gamma(δ^2/Vh, Vh/δ), β1, β2)
#
#function randlocus(M::AWGamma)
#    s = rand(M.sd)
#    g = rand(M.hd)
#    h = g - mean(M.hd) + M.β1/(1 + M.β2*s)
#    DipLocus(-s*h, -s)
#end
#
#function pdfh(dfe::AWGamma, h)
#    @unpack sd, hd, β1, β2 = dfe
#    d = mean(hd)
#    I, _ = quadgk(s->pdf(sd, s)*pdf(hd, h + d - β1/(1+β2*s)), 0, 1.)
#    return I
#end
