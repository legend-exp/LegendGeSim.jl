```@meta
EditURL = "simulate_realistic_waveforms_lit.jl"
```

# Simulate realistic waveforms

You can also download this tutorial as a
[Jupyter notebook](simulate_realistic_waveforms.ipynb) and a plain
[Julia source file](simulate_realistic_waveforms.jl).

Here we will use an example `pet` file `csv` format from `LegendTestData` corresponding to the Public Inverted Coax.

````@example simulate_realistic_waveforms
using LegendGeSim
using Plots
````

## Get inputs from `legend-test-data`

````@example simulate_realistic_waveforms
using LegendTestData
````

### Detector metadata

````@example simulate_realistic_waveforms
ldsim_path = joinpath(legend_test_data_path(), "data", "ldsim")

detector_name = "invcoax-metadata"
detector_metadata_filename = joinpath(ldsim_path, detector_name*".json");
nothing #hide
````

Alternatively, enter your own path to a real LEGEND detector JSON

```julia
detector_metadata_filename = "path/to/V04545A.json"
```

### PET input file

````@example simulate_realistic_waveforms
path_to_pet_file = joinpath(ldsim_path, "single-invcoax-th228-geant4.csv");
nothing #hide
````

## Settings

See manual on [Field Simulation](@ref) and [Ideal Pulse Simulation](@ref) for a detailed explanation of environment and simulation settings, as well as the noise model settings

````@example simulate_realistic_waveforms
environment_settings = Dict(
    "crystal_temperature_in_K" => 77,
    "medium" => "vacuum",
);
nothing #hide
````

simple settings for point charge simulation with dummy constant impurity

````@example simulate_realistic_waveforms
simulation_settings = Dict(
    "method" => "SSD",
    "cached_name" => "test",
);
nothing #hide
````

````@example simulate_realistic_waveforms
daq_settings = Dict(
    "preamp" => Dict(
        "type" => "generic",
        "t_decay_in_us" => 43, # from V04545A HADES data
        "t_rise_in_ns" => 100, # by eye
        "gain_ADC_eV" => 0.0138, # by eye from V04545A HADES data FEP @ 36100 ADC
        "offset_in_ADC" => 11900, # from V04545A HADES data mean() of baseline
        "noise_sigma_in_keV" => 2 # by eye
    ),
    "fadc" => Dict(
        "type" => "generic",
        "sampling_interval" => 16 # ns, from HADES data
    ),
    "trigger" => Dict(
        "type" => "trapezoidal",
        "window_lengths" => [250,250,250],
        "threshold" => 9 # keV
    ),
    "daq" => Dict(
        "type" => "generic",
        "nsamples" => 3748, # from HADES data
        "baseline_length" => 1770 # by eye from data
    )
);
nothing #hide
````

````@example simulate_realistic_waveforms
noise_model = Dict(
    "type" => "sim"
);
nothing #hide
````

## Simulate from scratch (`pet` -> `raw`)

````@example simulate_realistic_waveforms
raw_table = LegendGeSim.simulate_raw(detector_metadata_filename, path_to_pet_file, environment_settings, simulation_settings, daq_settings, noise_model; n_waveforms=10)
````

````@example simulate_realistic_waveforms
plot(raw_table.waveform)
````

The file contains double the amount of waveforms (20) compared to what we asked to simulate (`n_waveforms = 10`). That's because at the pulse simulation level, n+ contact pulses are kept as well. The first 10 entries are p+ contact waveforms.

````@example simulate_realistic_waveforms
length(raw_table)
````

````@example simulate_realistic_waveforms
plot(raw_table.waveform[1:10], legend=false, title="only p+ contact waveforms")
````

You can save simulated waveforms in a file that can be used later

````@example simulate_realistic_waveforms
using HDF5
using LegendHDF5IO
````

````@example simulate_realistic_waveforms
raw_name = "cache/test_100wfs_raw.hdf5"
h5open(raw_name, "w") do f
    LegendHDF5IO.writedata(f, "raw", raw_table[1:10])
end
````

## Simulate from `pss` table in code

Rather that simulating from scratch pet->raw, you may also input an already ready pss table with ideal pulses, and simulate only the pss->raw step, i.e. the DAQ chain simulation

````@example simulate_realistic_waveforms
pss_table, pss_truth = LegendGeSim.simulate_pulses(detector_metadata_filename, path_to_pet_file, environment_settings, simulation_settings, noise_model; n_waveforms=10);
nothing #hide
````

````@example simulate_realistic_waveforms
plot(pss_table.waveform[1:10])
````

````@example simulate_realistic_waveforms
raw_table1 = LegendGeSim.pss_to_raw(pss_table, pss_truth, daq_settings, noise_model)
````

````@example simulate_realistic_waveforms
plot(
    plot(pss_table.waveform[1:10]),
    plot(raw_table1.waveform[1:10]),
    size=(800,400)
)
````

Save the pss file for the next section

````@example simulate_realistic_waveforms
using LegendHDF5IO
using HDF5
````

````@example simulate_realistic_waveforms
pss_name = "cache/test_100wfs_pss.hdf5"
h5open(pss_name, "w") do f
    LegendHDF5IO.writedata(f, "pss/pss", pss_table[1:10])
    LegendHDF5IO.writedata(f, "pss/truth", pss_truth[1:10])
end
````

## Simulate from pre-saved `pss` hdf5 file

Rather than simulating form scratch pet->raw, you may input the name of a pre-saved pss file containing pss and pss truth information

````@example simulate_realistic_waveforms
raw_table2 = LegendGeSim.pss_to_raw(pss_name, daq_settings, noise_model)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

