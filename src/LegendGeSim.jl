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
using CurveFit
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
using LegendTextIO
# using LegendTextIO # Geant4CSVInput
using LinearAlgebra
using LsqFit
using MJDSigGen
# using ProgressMeter
using Parameters
using Polynomials
using PropDicts
using RadiationDetectorDSP # dev branch
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

# using Plots


T = Float32

const energy_unit = u"keV"
const ns_unit = u"ns"
const μs_unit = u"μs"

const germanium_ionization_energy = SolidStateDetectors.material_properties[:HPGe].E_ionisation # in eV

Random.seed!(123) # only for testing


include("sim_config.jl")
include("pss.jl")
include("detector.jl")

include("preamp.jl")
include("fadc.jl")

include("noise.jl")
include("baselines.jl")

include("elec_chain.jl")

include("trigger.jl")
include("daq.jl")

include("pet_to_stp.jl")
include("stp_to_pss.jl")
include("pss_to_raw.jl")
include("pet_to_raw.jl")

include("waveform_utils.jl")
include("temp_utils.jl")



end # module
