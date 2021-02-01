# This file is a part of LegendGeSim.jl, licensed under the MIT License (MIT).

__precompile__(true)

"""
    LegendGeSim

Template for Julia packages.
"""
module LegendGeSim

using LinearAlgebra
using Random
using Statistics

using ArgCheck
using ArraysOfArrays
using Clustering
using Distributions
using ElasticArrays
using EncodedArrays
using FillArrays
using IntervalSets
using JSON
using LegendDataTypes
using LegendHDF5IO
using RadiationDetectorSignals
using RadiationSpectra
using Random123
using RecipesBase
using Requires
using SolidStateDetectors
using StaticArrays
using StatsBase
using StructArrays
using Tables
using TypedTables
using Unitful


include("utils.jl")
include("g4_to_mcstp.jl")
include("mcstp_to_mcpss.jl")
include("mcpss_to_mcraw.jl")
# include("another_source_file.jl")
# ... more includes ...

end # module
