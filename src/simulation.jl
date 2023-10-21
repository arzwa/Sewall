
simulate!(args...; kwargs...) = simulate!(Random.default_rng(), args...; kwargs...)

function simulate!(rng::AbstractRNG, model, pop, n::Int; show_progress=true)
    prog = Progress(n; enabled=show_progress)
    for _=1:n
        generation!(rng, model, pop)
        next!(prog)
    end
    return pop
end

function simulate!(rng::AbstractRNG, model, pop, n, f::Function; 
        drop=0, thin=1, show_progress=true)
    x  = f(pop, model, 0)
    xs = Array{typeof(x)}(undef, n+1); xs[1] = x
    prog = Progress(n; enabled=show_progress)
    for i=1:n
        generation!(rng, model, pop)
        xs[i+1] = f(pop, model, i)
        next!(prog)
    end
    return pop, xs[drop+1:thin:end]
end

