"""
    Sewall

The goal of this module is to provide fast simulation methods for forward-time
population genetics simulation, specifically for subdivided populations and
biallelic finite-sites genetic architectures.
"""
module Sewall

using Printf
using QuadGK
using Random
using Reexport
using Parameters
using StatsBase
using Distributions
using LinearAlgebra
using ProgressMeter
@reexport using WrightDistribution
import WrightDistribution: expectedpq

include("architecture.jl")
include("deme.jl")
include("mainlandisland.jl")
include("metapopulation.jl")
include("simulation.jl")
include("effectivemigration.jl")
include("utils.jl")

vvcat(x) = vcat(x...)

export Locus, Architecture, Deme, FiniteIslands, MainlandIsland, FixedMainland
export generation!, initpop, vvcat, sasb, sehe, Ne2N, eqpdf, simulate!

end # module Sewall
