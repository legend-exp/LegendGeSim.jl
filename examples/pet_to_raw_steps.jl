using LegendGeSim
using HDF5
using LegendHDF5IO
filename(path) = splitext(basename(path))[1]

# filename(path) = splitext(basename(path))[1]
processed_dir = "output/"

## [WIP] 1. Load test Geant4 simulation of position-energy-timing information (PET)

## add LegendTestData module
# import Pkg
# Pkg.add(url="https://github.com/legend-exp/LegendTestData.jl")
# download test data
# Pkg.build("LegendTestData")

# using LegendTestData
# testdata_path = LegendTestData.legend_test_data_path()

testdata_path = "data"

## get path and name of input Geant4 CSV
# # position-energy-time information needed for pulse simulation
pet_input = "dual-invcoax-th228-geant4_small.csv"
# pet_input_fullpath = joinpath(testdata_path, "data/geant4", pet_input)
pet_input_fullpath = joinpath(testdata_path, pet_input)
det_meta = "public_ivc.json"
# det_meta_fullpath = joinpath(testdata_path, "data/json", det_meta)
det_meta_fullpath = joinpath(testdata_path, det_meta)


## simulation config 

# sim_config_filename = "configs/SSD_NoiseSim.json"
sim_config_filename = "configs/siggen_NoiseSim.json"


## ----- pet -> stp

##
# stp_table = LegendGeSim.pet_to_stp(pet_full_name, sim_config_filename)
stp_table = LegendGeSim.pet_to_stp(pet_input_fullpath, det_meta_fullpath, sim_config_filename)

## ------- intermediate step
# you can save it to use later, or skip saving and use mcstp_table in memory for next steps
stp_name = joinpath(processed_dir, "$(filename(pet_input_fullpath))_$(filename(det_meta_fullpath))_stp.h5")
HDF5.h5open(stp_name, "w") do f
    LegendHDF5IO.writedata(f, "stp", stp_table)
end

##

using Unitful
using Plots
histo_edep = histogram(sum.(ustrip.(stp_table.edep)))
png(histo_edep, "stp_edep.png")



# ----- stp -> pss

# pss_table, pss_truth = LegendGeSim.stp_to_pss(stp_table, det_meta_fullpath, sim_config_filename)
pss_table, pss_truth = LegendGeSim.stp_to_pss(stp_name, det_meta_fullpath, sim_config_filename)

## ------- intermediate step
# Saving table
pss_name = joinpath(processed_dir, "$(filename(pet_input_fullpath))_$(filename(det_meta_fullpath))_$(filename(sim_config_filename))_pss.h5")

h5open(pss_name, "w") do f
    LegendHDF5IO.writedata(f, "pss/pss", pss_table)
    LegendHDF5IO.writedata(f, "pss/truth", pss_truth)
end
@info "-> $(pss_name)"


##
using Plots
plot_pss = plot(pss_table.waveform[1:10])
# png(plot_pss, joinpath("plots", "$(out_filename_base)_pss.png"))

## ----- pss -> raw"

# raw_table = LegendGeSim.pss_to_raw(pss_table, pss_truth, det_meta_fullpath, sim_config_filename) 
raw_table = LegendGeSim.pss_to_raw(pss_name, det_meta_fullpath, sim_config_filename) 

##
using Plots
plot_raw = plot(raw_table.waveform[1:10], legend=false)#, linewidth=2)
# png(plot_raw, joinpath("plots", "$(out_filename_base)_raw.png"))

## Saving table
raw_name = joinpath(processed_dir, "$(filename(pet_input_fullpath))_$(filename(det_meta_fullpath))_$(filename(sim_config_filename))_raw.h5")
HDF5.h5open(raw_name, "w") do f
    LegendHDF5IO.writedata(f, "raw/raw", raw_table)
    LegendHDF5IO.writedata(f, "raw/truth", pss_truth)
end
@info "-> $raw_name"

## All steps together
pet_full_name = joinpath(pet_path, pet_name * ".hdf5")
raw_table = LegendGeSim.pet_to_raw(pet_full_name, sim_config_filename)