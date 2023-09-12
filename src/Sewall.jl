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
@reexport using WrightDistribution

include("architecture.jl")
include("deme.jl")
include("metapopulation.jl")

vvcat(x) = vcat(x...)

export Locus, Architecture, Deme, FiniteIslands
export generation!, initpop, vvcat, sasb, sehe, Ne2N

end # module Sewall
