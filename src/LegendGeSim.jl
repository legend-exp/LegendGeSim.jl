# This file is a part of LegendGeSim.jl, licensed under the MIT License (MIT).

__precompile__(true)

"""
    LegendGeSim

Template for Julia packages.
"""
module LegendGeSim

using ArgCheck
using ArraysOfArrays
using Clustering
using DataFrames # remove and use only Tables 
using Distributions
using DSP # Julia's Digital Signal Processor (DSP) Package
using ElasticArrays
using EncodedArrays
using FillArrays
using HDF5 # not in Project.toml and not supposed to be (LegendHDF5IO)
using IntervalSets
using JSON
using LegendDataTypes
using LegendHDF5IO
# using LegendTextIO # Geant4CSVInput
using LinearAlgebra
# using RadiationDetectorDSP # By Oliver Schulz
# using ProgressMeter
using Parameters
using PropDicts
using RadiationDetectorSignals
using RadiationSpectra
using Random
using Random123
using RecipesBase
using Requires
using SolidStateDetectors
using StaticArrays
using Statistics
using StatsBase
using StructArrays
using Tables
using TypedTables
using Unitful



include("utils.jl")
include("elec.jl")
include("g4_to_mcstp.jl")
include("mcstp_to_mcpss.jl")
include("mcpss_to_mcraw.jl")
include("g4_to_mcraw.jl")
# include("another_source_file.jl")
# ... more includes ...

end # module
