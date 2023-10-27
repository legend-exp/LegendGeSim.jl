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
    "cached_name" => "test",
);

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

noise_model = Dict(
    "type" => "sim"
);

raw_table = LegendGeSim.simulate_raw(detector_metadata_filename, path_to_pet_file, environment_settings, simulation_settings, daq_settings, noise_model; n_waveforms=10)

plot(raw_table.waveform)

length(raw_table)

plot(raw_table.waveform[1:10], legend=false, title="only p+ contact waveforms")

using HDF5
using LegendHDF5IO

raw_name = "cache/test_100wfs_raw.hdf5"
h5open(raw_name, "w") do f
    LegendHDF5IO.writedata(f, "raw", raw_table[1:10])
end

pss_table, pss_truth = LegendGeSim.simulate_pulses(detector_metadata_filename, path_to_pet_file, environment_settings, simulation_settings, noise_model; n_waveforms=10);

plot(pss_table.waveform[1:10])

raw_table1 = LegendGeSim.pss_to_raw(pss_table, pss_truth, daq_settings, noise_model)

plot(
    plot(pss_table.waveform[1:10]),
    plot(raw_table1.waveform[1:10]),
    size=(800,400)
)

using LegendHDF5IO
using HDF5

pss_name = "cache/test_100wfs_pss.hdf5"
h5open(pss_name, "w") do f
    LegendHDF5IO.writedata(f, "pss/pss", pss_table[1:10])
    LegendHDF5IO.writedata(f, "pss/truth", pss_truth[1:10])
end

raw_table2 = LegendGeSim.pss_to_raw(pss_name, daq_settings, noise_model)
