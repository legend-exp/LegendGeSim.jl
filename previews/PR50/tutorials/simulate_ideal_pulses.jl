using LegendGeSim
using Plots

using LegendTestData

ldsim_path = joinpath(legend_test_data_path(), "data", "ldsim")

detector_name = "invcoax-metadata"
detector_metadata_filename = joinpath(ldsim_path, detector_name*".json");

path_to_pet_file = joinpath(ldsim_path, "single-invcoax-th228-geant4.csv");

environment_settings = Dict(
    "crystal_temperature_in_K" => 77,
    "medium" => "vacuum",
);

simulation_settings = Dict(
    "method" => "SSD",
    "cached_name" => "", # a non-empty string will cache the simulation results
);

noise_model = Dict(
    "type" => "sim"
);

# Simulate ideal pulses

pss_table, pss_truth = LegendGeSim.simulate_pulses(detector_metadata_filename, path_to_pet_file, environment_settings, simulation_settings; n_waveforms=5);

plot(pss_table.waveform, legend=false, linewidth=1.5)

length(pss_table)

plot(pss_table.waveform[1:5], legend=false, title="only p+ contact pulses")

savefig("ideal_pulses.pdf") # hide

using HDF5
using LegendHDF5IO

pss_name = "cache/test_100wfs_pss.hdf5"
h5open(pss_name, "w") do f
    LegendHDF5IO.writedata(f, "pss/pss", pss_table[1:5])
    LegendHDF5IO.writedata(f, "pss/truth", pss_truth[1:5])
end
