using LegendGeSim
using HDF5
using LegendHDF5IO

## [WIP] 1. Load test Geant4 simulation of position-energy-timing information (PET)

## add LegendTestData module
# import Pkg
# Pkg.add(url="https://github.com/legend-exp/LegendTestData.jl")
# download test data
# Pkg.build("LegendTestData")

## get path and name of input Geant4 CSV
# using LegendTestData
# # position-energy-time information needed for pulse simulation
# pet_name = "dual-invcoax-th228-geant4.csv"
# pet_path = LegendTestData.legend_test_data_path()
# processed_dir = "cache/"

## HADES
pet_name = "raw-IC160A-Th228-uncollimated-top-run0002-source_holder-bi-hdf5-01-test"
pet_path = "data/"
processed_dir = "cache/"


## simulation config 

# sim_config_filename = "data/V02160A_siggen_NoiseData.json"
# sim_config_filename = "data/V02160A_siggen_NoiseSim.json"
sim_config_filename = "data/V02160A_SSD_NoiseSim.json"
# sim_config_filename = "data/V02160A_SSD_NoiseData.json"

filename(path) = splitext(basename(path))[1]
out_filename_base = filename(sim_config_filename)

## ----- pet -> stp
pet_full_name = joinpath(pet_path, pet_name)

## test 
pet_file = LegendGeSim.read_pet(pet_full_name)

##
stp_table = LegendGeSim.pet_to_stp(pet_full_name, sim_config_filename)

## ------- intermediate steps
# you can save it to use later, or skip saving and use mcstp_table in memory for next steps
stp_name = joinpath(processed_dir, out_filename_base*"_stp.h5")
HDF5.h5open(stp_name, "w") do f
    LegendHDF5IO.writedata(f, "stp", stp_table)
end

## Read stp (if the previous step was skipped, otherwise skip this step)
stp_name = joinpath(processed_dir, out_filename_base*"_stp.h5")
stp_table = HDF5.h5open(stp_name, "r") do input
    LegendHDF5IO.readdata(input, "stp")
end

##

using Unitful
using Plots
histo_edep = histogram(sum.(ustrip.(stp_table.edep)))
png(histo_edep, "stp_edep.png")

# ------- intermediate steps


# ----- stp -> pss

pss_table, pss_truth = LegendGeSim.stp_to_pss(stp_table, sim_config_filename)

## ------- intermediate steps
# Saving tablec
# noise_name = occursin("Data", sim_config_filename) ? "NoiseData" : "NoiseSim"
# out_filename = joinpath(processed_dir, "$(mc_name)_$(det_name)_$(sim_config.sim_method)_$(noise_name)_mcpss.h5")
out_filename = joinpath(processed_dir, out_filename_base*"_pss.h5")
h5open(out_filename, "w") do f
    LegendHDF5IO.writedata(f, "pss/pss", pss_table)
    LegendHDF5IO.writedata(f, "pss/truth", pss_truth)
end
@info "-> $(out_filename)"

## Read pss (if the previous step was skipped, otherwise skip this step)

pss_name = joinpath(processed_dir, out_filename_base*"_pss.h5")

pss_h5 = HDF5.h5open(pss_name, "r")
pss_table = LegendHDF5IO.readdata(mcpss_h5, "pss/pss")
pss_truth = LegendHDF5IO.readdata(mcpss_h5, "pss/truth")
HDF5.close(pss_h5)

# ------- intermediate steps


##
using Plots
plot_pss = plot(pss_table.waveform[1:10])
png(plot_pss, joinpath("plots", "$(out_filename_base)_pss.png"))

## ----- pss -> raw"

raw_table = LegendGeSim.pss_to_raw(pss_table, pss_truth, sim_config_filename) 

##
using Plots
plot_raw = plot(raw_table.waveform[1:10], legend=false)#, linewidth=2)
png(plot_raw, joinpath("plots", "$(out_filename_base)_raw.png"))

## Saving table
out_filename = joinpath(processed_dir, out_filename_base*"_raw.h5")
HDF5.h5open(out_filename, "w") do f
    LegendHDF5IO.writedata(f, "raw/raw", raw_table)
    LegendHDF5IO.writedata(f, "raw/truth", pss_truth)
end
@info "-> $out_filename"

## All steps together
pet_full_name = joinpath(pet_path, pet_name * ".hdf5")
raw_table = LegendGeSim.pet_to_raw(pet_full_name, sim_config_filename)