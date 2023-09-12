using Sewall
using Random; rng = Random.default_rng()
using Plots, PlotThemes; theme(:hokusai)


# Reduction in diversity due to a sweep
L = 50
N = 500
k = 1
s = 0.01
u = s/20
r = 0.02
x = [Locus(0.0, 0.0, 0.0) for _=1:(L+1)]
j = L÷2 +1 
A = Architecture(x, fill(u, L+1), fill(r, L))
A.u[j] = 0.0
D = Deme(N=N, k=k, A=A)
P = initpop(rng, D, zeros(L+1)) 
ps = map(1:11000) do _
    generation!(rng, D, P)
    sum(P.haploids, dims=1) ./ N
end
while true
    A.loci[j] = Locus(0.1, 0.0, 0.0)
    A.u[j] = 0.0
    P.haploids[1,L÷2+1] = true
    ys = map(1:500) do _
        generation!(rng, D, P)
        sum(P.haploids, dims=1) ./ N
    end
    if ys[end][j] == 1.0 
        push!(ps, ys...)
        break
    end
end
xs = vcat(ps...)

pp = xs .* (1 .- xs)
plot(pp[:,1])
plot!(pp[:,10])
plot!(pp[:,25])

# BGS
L  = 80
N  = 1000
k  = 1
Ns = 10.0
s  = Ns * (1/N + (1/(2N*k)))
js = 30:50
ud = 5/(N*length(js))
u  = ud
r  = s/500
x  = [Locus(0.0, 0.0, 0.0) for _=1:L]
A  = Architecture(x, fill(u, L+1), fill(r, L))
for j=js; A.loci[j] = Locus(-s, 0.0, 0.0); end
A.u[js] .= ud
D = Deme(N=N, k=k, A=A)
P = initpop(rng, D, zeros(L)) 
ps = map(1:110000) do n
    n % 100 == 0 && (@info n)
    generation!(rng, D, P)
    sum(P.haploids, dims=1) ./ N
end |> vvcat

plot(map(var, eachcol(ps[1001:100:end,:])))

ds = Sewall.eqpdf(D)
stephist(ps[:,30], bins=0:0.01:1, norm=:probability)
xs, ys = sfs(ds[30], 0:0.01:1)
plot!(xs, reverse(ys ./ sum(ys)))


# Metapopulation
d = 10
model = FiniteIslands(demes=[D for i=1:d], 
                      mhap=rand(d) ./ 5, 
                      mdip=rand(d) ./ 5)

pops = initpop(rng, model, [rand(L) for i=1:d])
generation!(rng, model, pops)

for i=1:1000
    i % 100 == 0 && (@info i)
    generation!(rng, model, pops)
end


mat = @view P.diploids[1,:]
pat = @view P.diploids[2,:]
tgt = @view P.haploids[1,:]

Sewall.meiosis!(rng, tgt, A, mat, pat) 

P.haploids[1,1] = false

#----
s = 0.05
u = s/5
A = Architecture([Locus(-s, 0.0, 0.0)], [u], Float64[])
D = Deme(N=220, k=2, A=A)
ps = map(1:200) do _
    P = initpop(rng, D, rand(length(A))) 
    for _=1:1000
        generation!(rng, D, P)
    end 
    sum(P.haploids) ./ D.N
end

f, = Sewall.eqpdf(D)
bs = 0:0.02:1
stephist(ps, bins=bs, norm=:probability)
ys = map(p->pdf(f, 1-p), bs)
plot!(bs, ys ./sum(ys))
vline!([u/s])

# BGS proof of principle
# ======================
N = 3000
k = 1
Ne = 1/(1/N + 1/(2N*k))
s  = 100/Ne
u = 0.4*s
r = 0.001*s
u0 = s/100
# prediction based on Nordborg 1997
NeBGS(Ne, u, s, r) = (1 - (u/s)/((1+(r/s))^2))*Ne
Ne_ = NeBGS(Ne, u, s, r)
d1 = Wright(Ne,  u0, u0, 0., 0.)
d2 = Wright(Ne_, u0, u0, 0., 0.)
bins = 0:0.02:1
y1 = map(p->pdf(d1, p), bins)
y2 = map(p->pdf(d2, p), bins)
plot(bins, y1)
plot!(bins, y2)

A = Architecture([
    Locus(-s, 0.0, 0.0, u, 0.0), 
    Locus(0.0, 0.0, 0.0, u0, u0)], [r]) 
D = Deme(N=N, k=k, A=A)
P = initpop(rng, D, zeros(2))
    
ngen = 500ceil(Int, Ne)
ps = map(1:ngen) do i
    i % 5000 == 0 && (@info "$i/$ngen")
    generation!(rng, D, P)
    sum(P.haploids, dims=1) / D.N
end |> vvcat

stephist(ps[5ceil(Int,Ne):100:end,2], bins=bins, norm=true, fill=true,
    fillalpha=0.3, label="simulation", xlabel="\$p\$", ylabel="\$\\phi\$")
plot!(bins, y1, label="\$N_e\$")
plot!(bins, y2, label="\$N_{e,BGS}\$", size=(500,300))

