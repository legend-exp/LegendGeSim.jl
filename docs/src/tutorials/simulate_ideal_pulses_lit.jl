# # Simulate Ideal Pulses

#md # You can also download this tutorial as a
#md # [Jupyter notebook](simulate_ideal_pulses.ipynb) and a plain
#md # [Julia source file](simulate_ideal_pulses.jl).

# Here we will use an example `pet` file `csv` format from `LegendTestData` corresponding to the Public Inverted Coax.

using LegendGeSim
using Plots

# ## Get inputs from `legend-test-data`

using LegendTestData

# ### Detector metadata

ldsim_path = joinpath(legend_test_data_path(), "data", "ldsim")

detector_name = "invcoax-metadata"
detector_metadata_filename = joinpath(ldsim_path, detector_name*".json");

# Alternatively, enter your own path to a real LEGEND detector JSON

# ```julia
# detector_metadata_filename = "path/to/V04545A.json"
# ```

# ### PET input file

path_to_pet_file = joinpath(ldsim_path, "single-invcoax-th228-geant4.csv");

# ## Settings

# See manual on [Field Simulation](@ref) and [Ideal Pulse Simulation](@ref) for a detailed explanation of environment and simulation settings, as well as the noise model settings

environment_settings = Dict(
    "crystal_temperature_in_K" => 77, 
    "medium" => "vacuum", 
);

# Simple settings for point charge simulation with dummy constant impurity

simulation_settings = Dict(
    "method" => "SSD",
    "cached_name" => "test",
);

#

noise_model = Dict(
    "type" => "sim"
);

## Simulate ideal pulses

# Provide an optional argument `n_waveforms` to define the number of pulses to be simulated. Default: all pulses based on the input file

pss_table, pss_truth = LegendGeSim.simulate_pulses(detector_metadata_filename, path_to_pet_file, environment_settings, simulation_settings; n_waveforms=5);

# 

plot(pss_table.waveform, legend=false, linewidth=1.5)

# The "negative" waveforms are the ones coming from the n+ contact, so they are symmetric. In fact, the length of our table is twice the amount of waveforms we simulated.

length(pss_table)

# 

plot(pss_table.waveform[1:5], legend=false, title="only p+ contact pulses")

#jl savefig("ideal_pulses.pdf") # hide
#md savefig("ideal_pulses.svg"); nothing # hide
#md # ![ideal_pulses](ideal_pulses.svg)

# ## Save output to hdf5 file

# You can save simulated pulses in a file that can be used later.

using HDF5
using LegendHDF5IO

#

pss_name = "cache/test_100wfs_pss.hdf5"
h5open(pss_name, "w") do f
    LegendHDF5IO.writedata(f, "pss/pss", pss_table[1:5])
    LegendHDF5IO.writedata(f, "pss/truth", pss_truth[1:5])
end

