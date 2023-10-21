using Sewall
using Test, Random; rng = Random.default_rng()

@testset "Mutation-drift balance" begin
    A = Architecture([Locus(0.0, 0.0, 0.0, 0.1, 0.1)], Float64[])
    D = Deme(N=320, k=5, A=A)
    ps = map(1:200) do _
        P = initpop(rng, D, rand(length(A))) 
        for _=1:200
            generation!(rng, D, P)
        end 
        sum(P.haploids) ./ D.N
    end
    # standard error on the mean allele frequency should be 
    # √((p(1-p)/Ne)/nrep)
    se = √(0.5 / D.Ne / length(ps))
    @test abs(mean(ps) - 0.5) < 2se
    # what would be a statistical test for the variance?
end

# Test mutation-selectioni-drift balance
@testset "Mutation-selection-drift balance" begin
    s = 0.05
    u = s/10
    A = Architecture([Locus(-s, 0.0, 0.0, u, u)], Float64[])
    D = Deme(N=220, k=2, A=A)
    ps = map(1:200) do _
        P = initpop(rng, D, rand(length(A))) 
        for _=1:1000
            generation!(rng, D, P)
        end 
        sum(P.haploids) ./ D.N
    end
    # expected value is u/s=0.1, variance should be ≈ p(1-p)/Ne... but the
    # deviation is quite substantial...
    # standard error on the mean allele frequency should be √((p(1-p)/Ne)/nrep)
    d,  = Sewall.eqpdf(D)
    se = √(var(d) / length(ps))
    @test abs(mean(ps) - u/s) < 2se
end

# slight upward bias? or is this expected due to drift?
#@testset "Quasi-deterministic mutation-selection balance" begin
#    s = 0.05
#    u = s/10
#    A = Architecture([Locus(-s, 0.0, 0.0, u, 0.0)], Float64[])
#    D = Deme(N=30000, k=1, A=A)
#    P = initpop(rng, D, rand(length(A))) 
#    ps = Float64[]
#    for i=1:60000
#        generation!(rng, D, P)
#        if i%1000 == 0 
#            @info i; push!(ps, sum(P.haploids) / D.N)
#        end
#    end 
#    @test abs(mean(ps) - u/s) < 0.01
#end

# Check the decay of LD.
@testset "Two-locus LD" begin
    ld(X) = cov(X)[2,1]
    for r=[0.5, 0.1, 0.05, 0.01]
        A = Architecture([Locus(0.0, 0.0, 0.0, 0.0, 0.0) for _=1:2], [r])
        D = Deme(N=100000, k=1, A=A)
        P = initpop(rng, D, zeros(2)) 
        for i=1:(D.N÷2)
            P.haploids[i,:] .= [true,true]
            P.haploids[D.N÷2 + i,:] .= [false,false]
        end
        ds = map(0:10) do _
            d = ld(P.haploids)
            generation!(rng, D, P); d
        end 
        eds = map(t->0.25*(1-r)^t, 0:10)
        @test all(abs.(ds .- eds) .< 0.01) 
    end
end

@testset "Two-locus LD in L locus architecture" begin
    L = 10
    ld(X) = cov(X)[1,L]
    for r=[0.5, 0.1, 0.05, 0.01]
        A = Architecture(Locus(0.0, 0.0, 0.0, 0.0, 0.0), L, r)
        D = Deme(N=100000, k=1, A=A)
        P = initpop(rng, D, zeros(L)) 
        for i=1:(D.N÷2)
            P.haploids[i,:] .= ones(Bool, L)
            P.haploids[D.N÷2 + i,:] .= zeros(Bool, L)
        end
        ds = map(0:10) do _
            d = ld(P.haploids)
            generation!(rng, D, P); d
        end 
        eds = map(t->0.25*(1-D.A.R[1,L])^t, 0:10)
        @test all(abs.(ds .- eds) .< 0.01) 
        # ds = round.(ds, digits=3)
        # @show collect(zip(ds, eds))
    end
end

@testset "Mainland-island model" begin
    s = 0.05
    u = s/200
    A = Architecture([Locus(-s, 0.0, 0.0, u, u)], Float64[])
    D = Deme(N=220, k=5, A=A)
    y = FixedMainland([true], ([true],[true]))
    M = MainlandIsland(D, s*0.2, 0.0, y)
    ps = map(1:10) do _
        P = initpop(rng, M, rand(length(A))) 
        map(1:20000) do _
            generation!(rng, M, P)
        end 
        sum(P.haploids) ./ D.N
    end
    d, = eqpdf(M)
    se = √(var(d) / length(ps))
    @test abs(mean(ps) - (1-mean(d))) < 2se
end

