function windowsmooth(x::Vector{T}, y::Vector{T}, window, overlap, dropna=true) where T
    imin = 1
    xs = T[]
    ys = T[]
    @info "" collect(0:overlap:x[end]-window)
    for x0=0:overlap:x[end]-window
        idx = [i for i=imin:length(x) if x0 <= x[i] < x0 + window]
        push!(xs, x0 + window/2)
        if length(idx) > 0 
            push!(ys, mean(y[idx]))
            imin = idx[1]
        else
            push!(ys, NaN)
        end
    end
    return (xs, ys)
end

