using Sewall, Plots

PP = map(enumerate([0., 0.5, 1.0])) do (j,h)
    P = plot()
    xmx = 1.25
    map(enumerate([0.01,0.25,0.5,0.75,1.,1.25,1.5])) do (i,Ls)
        s  = 0.01
        sa = -s*h
        sb = -s -2sa
        x, y, x_, y_, = Sewall.findroots_ms(sa, sb, Ls/s)
        color = i == 1 ? :lightgray : i-1
        plot!(x,y,   color=color, lw=3, alpha=0.7, label="\$Ls=$(@sprintf "%.2f" Ls)\$")
        plot!(x_,y_, color=color, alpha=0.7, label="", ls=:dot)
        plot!([last(x), last(x), xmx], [last(y), 0, 0], color=color, lw=3, alpha=0.7, label="")
    end
    lab = ["A", "B", "C"][j]
    P1 = plot(P, title="($lab) \$s=0.01, h=$(@sprintf "%.1f" h)\$", 
              xlabel="\$m/s\$", xlim=(0,xmx), bottom_margin=5Plots.mm,
              ylabel="\$\\tilde{p}\$", legend=h == 0.5 ? :bottomright : false)
end
plot(PP...)

