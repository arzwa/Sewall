using Literate
Literate.markdown(
    joinpath(@__DIR__, "README.jl"),
    joinpath(@__DIR__, "../"),
    documenter=false, execute=true)
