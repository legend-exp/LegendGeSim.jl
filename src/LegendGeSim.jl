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
using DelimitedFiles
using Distributions
using DSP # Julia's Digital Signal Processor (DSP) Package
using ElasticArrays
using EncodedArrays
using FillArrays
using HDF5 # later not supposed to be here, use LegendHDF5IO functions
using IntervalSets
using JSON
using LegendDataTypes
using LegendHDF5IO
# using LegendTextIO # Geant4CSVInput
using LinearAlgebra
using MJDSigGen
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
#
using CurveFit
using Polynomials
using LsqFit

T = Float32

const energy_unit = u"keV"
const ns_unit = u"ns"
const μs_unit = u"μs"

# germanium_ionization_energy = T(2.95)u"eV"
const germanium_ionization_energy = SolidStateDetectors.material_properties[:HPGe].E_ionisation # already in eV

Random.seed!(123) # only for testing

include("filters.jl")

include("sim_config.jl")
include("pss.jl")

include("detector.jl")

include("noise.jl")

include("elec_chain.jl")
include("daq.jl")

include("g4_to_mcstp.jl")
include("mcstp_to_mcpss.jl")
include("mcpss_to_mcraw.jl")

include("g4_to_mcraw.jl")

include("baselines.jl")

# move to LegendHDF5IO?
include("io.jl")

end # module
