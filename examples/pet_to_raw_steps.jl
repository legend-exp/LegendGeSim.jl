using LegendGeSim
using HDF5
using LegendHDF5IO
filename(path) = splitext(basename(path))[1]

processed_dir = "output/"

## [DO ONCE] Set up LegendTestData to get Geant4 simulation
# that serves as input position-energy-timing information (PET) for PSS

# import Pkg
# Pkg.add(url="https://github.com/legend-exp/LegendTestData.jl")
# download test data
# Pkg.build("LegendTestData")

# using LegendTestData
# testdata_path = LegendTestData.legend_test_data_path()
testdata_path = "data"

## select files of interest

## get path and name of input Geant4 CSV
# # position-energy-time information needed for pulse simulation
pet_input = "dual-invcoax-th228-geant4_small.csv"

# pet_input_fullpath = joinpath(testdata_path, "data/geant4", pet_input)
pet_input_fullpath = joinpath(testdata_path, pet_input)
det_meta = "public_ivc.json"
det_meta_fullpath = joinpath(testdata_path, det_meta)

## test

# pet_table = LegendGeSim.read_pet(pet_input_fullpath)
# histo_edep = histogram(ustrip.(pet_table.edep))

## simulation config 

# sim_config_filename = "configs/SSD_NoiseSim.json"
sim_config_filename = "configs/siggen_NoiseSim.json"


## ----- pet -> stp

stp_table = LegendGeSim.pet_to_stp(pet_input_fullpath, det_meta_fullpath)

## ------- intermediate steps
# you can save it to use later, or skip saving and use mcstp_table in memory for next steps
# stp_name = joinpath(processed_dir, out_filename_base*"_stp.h5")
# propdict(sim_config_filename).pss.cached_name
stp_name = joinpath(processed_dir, "$(filename(pet_input_fullpath))_$(filename(det_meta_fullpath))_stp.h5")
HDF5.h5open(stp_name, "w") do f
    LegendHDF5IO.writedata(f, "stp", stp_table)
end

# ## Read stp (if the previous step was skipped, otherwise skip this step)
# stp_name = joinpath(processed_dir, "$(filename(pet_input_fullpath))_$(filename(det_meta_fullpath))_stp.h5")
# stp_table = HDF5.h5open(stp_name, "r") do input
#     LegendHDF5IO.readdata(input, "stp")
# end

##

using Unitful
using Plots
histo_edep = histogram(sum.(ustrip.(stp_table.edep)))
# plot_name = joinpath("plots", "$(filename(pet_input_fullpath))_$(filename(det_meta_fullpath))_stp.png")
# png(histo_edep, plot_name)

# ------- intermediate steps


# ----- stp -> pss

pss_table, pss_truth = LegendGeSim.stp_to_pss(stp_table, det_meta_fullpath, sim_config_filename)
# pss_table, pss_truth = LegendGeSim.stp_to_pss(stp_name, det_meta_fullpath, sim_config_filename)

## ------- intermediate steps
# Saving tablec
# noise_name = occursin("Data", sim_config_filename) ? "NoiseData" : "NoiseSim"
# out_filename = joinpath(processed_dir, "$(mc_name)_$(det_name)_$(sim_config.sim_method)_$(noise_name)_mcpss.h5")
# out_filename = joinpath(processed_dir, out_filename_base*"_pss.h5")
pss_name = joinpath(processed_dir, "$(filename(pet_input_fullpath))_$(filename(det_meta_fullpath))_$(filename(sim_config_filename))_pss.h5")

h5open(pss_name, "w") do f
    LegendHDF5IO.writedata(f, "pss/pss", pss_table)
    LegendHDF5IO.writedata(f, "pss/truth", pss_truth)
end
@info "-> $(pss_name)"


# ------- intermediate steps


##
using Plots
plot_pss = plot(pss_table.waveform[1:10])
pss_plot = joinpath("plots", "$(filename(pet_input_fullpath))_$(filename(det_meta_fullpath))_$(filename(sim_config_filename))_pss.png")
png(plot_pss, pss_plot)

## ----- pss -> raw"

raw_table = LegendGeSim.pss_to_raw(pss_table, pss_truth, det_meta_fullpath, sim_config_filename) 
# raw_table = LegendGeSim.pss_to_raw(pss_name, det_meta_fullpath, sim_config_filename) 

##
using Plots
plot_raw = plot(raw_table.waveform[1:10], legend=false)#, linewidth=2)
raw_plot = joinpath("plots", "$(filename(pet_input_fullpath))_$(filename(det_meta_fullpath))_$(filename(sim_config_filename))_raw.png")
png(plot_raw, raw_plot)

## Saving table
raw_name = joinpath(processed_dir, "$(filename(pet_input_fullpath))_$(filename(det_meta_fullpath))_$(filename(sim_config_filename))_raw.h5")
HDF5.h5open(raw_name, "w") do f
    LegendHDF5IO.writedata(f, "raw/raw", raw_table)
    LegendHDF5IO.writedata(f, "raw/truth", pss_truth)
end
@info "-> $raw_name"

## All steps together
raw_table = LegendGeSim.pet_to_raw(pet_input_fullpath, det_meta_fullpath, sim_config_filename)
