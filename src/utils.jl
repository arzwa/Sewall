function windowsmooth(x, y, window, overlap, dropna=true)
    imin = 1
    res = map(0:overlap:x[end]) do x0
        i1 = findfirst(xj->xj > x0,        x[imin:end])
        i2 = findfirst(xj->xj > x0+window, x[imin:end])
        ys = isnothing(i2) ? 
            y[(imin + i1 - 1):end] :
            y[(imin + i1 - 1):(imin + i2 - 1)]
        ys = dropna ? filter(!isnan, ys) : ys
        x̄ = x0 + window/2
        imin += i1
        x̄, mean(ys)
    end
    first.(res), last.(res)
end

