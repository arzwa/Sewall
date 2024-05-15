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
using Roots
using ForwardDiff
using Reexport
using Parameters
using StatsBase
using Distributions
using LinearAlgebra
using ProgressMeter
@reexport using WrightDistribution
import WrightDistribution: expectedpq

include("architecture.jl")
include("recombination.jl")
include("deme.jl")
include("mainlandisland.jl")
include("metapopulation.jl")
include("simulation.jl")
include("effectivemigration.jl")
include("utils.jl")
include("qtrait.jl")
include("deterministic.jl")
include("dfe.jl")

vvcat(x) = vcat(x...)

export Locus, Architecture, Deme, FiniteIslands, MainlandIsland, FixedMainland, HWLEMainland
export generation!, initpop, vvcat, sasb, sehe, Ne2N, eqpdf, simulate!, equilibrium
export humanmap, humanmap2, flymap, flymap2

end # module Sewall
