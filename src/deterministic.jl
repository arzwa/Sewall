# Deterministic bifurcation analysis for the multilocus model
function tosolve(θ)
	@unpack sa, sb, m, L = θ
	g(p) = exp(2L*p*(sa + sb*(1-p)))
	f(p) = -m*g(p)*p - p*(1-p)*(sa + sb*(1-p))
end

function roots_stab(sa, sb, m, L)
	f  = tosolve((sa=sa, sb=sb, m=m, L=L))
	zs = find_zeros(f, 0.,  1.)
	filter!(!iszero, zs)
	ds = ForwardDiff.derivative.(f, zs)
	zs, ds
end

"""
    findroots_ms(sa, sb, L; stepsize, tol, maxit)

Find the stable and unstable (if any) internal fixed points for a deterministic
multilocus model with `L` equal effect loci. This finds the roots at intervals,
decreasing the step adaptively until they collide.
"""
function findroots_ms(sa, sb, L; ms0=0., stepsize=0.005, tol=1e-2, maxit=10000)
    s = -(sb + 2sa)  # this is the s scale
	ms = ms0
	zs, ds = roots_stab(sa, sb, ms*s, L)
	sol = [(ms, z, d) for (z, d) in zip(zs,ds)]
	i = 0
	while i < maxit && length(zs) >= 1 
		i += 1
		(length(zs) == 2 && abs(zs[1]-zs[2]) < tol) && break
		ms += stepsize
		zs_, ds = roots_stab(sa, sb, ms*s, L)
		if length(zs) == 2 && length(zs_) < 2 
			ms -= stepsize
			stepsize /= 3
		elseif length(zs) == 1 && length(zs_) == 1
			push!(sol, (ms, zs_[1], ds[1]))
            zs = zs_
		elseif length(zs_) == 2
			push!(sol, (ms, zs_[1], ds[1]))
			push!(sol, (ms, zs_[2], ds[2]))
            zs = zs_
		end
	end
	x = first.(sol)
	y = getindex.(sol, 2)
	z = last.(sol)
	stab = z .< 0
	unstab = z .> 0
	x[stab], y[stab], x[unstab], y[unstab] 
end

# find the critical m/s (either where two fixed points collide, or where the
# stable equilibrium disappears).
function critical_ms(sa, sb, L; mmax=1., tol=1e-2)
    s = -(sb + 2sa)  # this is the s scale 
    lb = 0.
    # find an upper bound
    ub = mmax
    zs, ds = roots_stab(sa, sb, ub*s, L)
    while length(zs) > 0
        lb = ub
        ub *= 2
        zs, ds = roots_stab(sa, sb, ub*s, L)
    end
    # bisect the interval that we found until we reach the desired tolerance
    ms = (ub + lb) / 2
    zs, ds = roots_stab(sa, sb, ms*s, L)
    while true
        if length(zs) > 0 
            lb = ms
        else
            ub = ms
        end
        ms′ = (ub + lb) / 2
        abs(ms′ - ms) < tol && break
        ms = ms′
        zs, ds = roots_stab(sa, sb, ms*s, L)
    end
    # the upper bound
    zs, ds = roots_stab(sa, sb, lb*s, L)
    return ms, zs, ds
end

