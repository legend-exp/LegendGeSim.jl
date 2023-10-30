# Realistic Waveform Simulation

As inputs to realistic waveform simulation you should provide

1. LEGEND detector metadata

2. PET input file

3. Environment settings

4. Simulation settings

5. DAQ settings

6. Noise model

## Example

```julia
detector_metadata_filename = "path/to/legend-detectors/germanium/diodes/V04545A.json"

pet_input_file = "/lfs/l1/legend/detector_char/enr/hades/simulations/legend-g4simple-simulation/simulations/V04545A/am_HS1/top_46r_4z/hdf5/sim-V04545A-am_HS1-top-46r-4z-01.hdf5"

environment_settings = Dict(
    "crystal_temperature_in_K" => 77, 
    "medium" => "vacuum", 
    "dl" => "vendor" # optional, default 0
    "operating_voltage_in_V" => 5000, # optional, default recV from metadata
)

simulation_settings = Dict(
    "method" => "SSD",
    "cached_name" => "test",
    "crystal_metadata_path" => "path/to/legend-detectors/germanium/crystals",
    "time_step" => 2, #ns
    "diffusion" => true,
    "self-repulsion" => true,
    "num_carriers" => 10
)

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
)

noise_model = Dict(
    "type" => "sim"
)

raw_table = LegendGeSim.simulate_raw(detector_metadata_filename, pet_input_file, environment_settings, simulation_settings, daq_settings, noise_model; n_waveforms=100)
```

The optional variable `n_waveforms` defines how many waveforms will be simulted based on the given `pet` file. Not providing the variable will default `0` meaning, simulate all waveforms in the `pet` input file.

The output `raw_table` is a Julia `Table`, currently mimicking HADES raw data format. It contains the following fields:
- `baseline`: baseline offset of each waveform in units of ADC
- `channel`: currently, ID of the contact (1 for p+ contact, 2 for n+ contact)
- `energy`: online energy estimated from the trap filter trigger. Currently in units of keV, data raw seems to be in ADC (WIP)
- `ievt`: event number
- `numtraces`: number of triggered detectors (always 1 for HADES)
- `packet_id`: means to packet losses (always 0 in simulation)
- `timestamp`: currently frist MC truth hit time of each event
- `tracelist`: lists of ADCs that triggered, 1 for HADES all the time
- `waveform`: waveform in format of `RDWaweform` from `RadiationDetectorSignals.jl`. Contains fields `time` and `signal`
- `wf_max`: maximum of each waveform
- `wf_std`: standard deviation of each waveform (full waveform not just baseline)

This format is digestible by pygama `build_dsp()` with HADES configuration JSON file.

## 1. LEGEND detector metadata

To obtain JSONs for all LEGEND detectors you can clone the repository `legend-detectors`. The JSONs can be found under `germanium/diodes`. You may also create your custom detector JSON file as long as it follows the LEGEND detector metadata format (see README in `diodes/`).

## 2. PET input file

See [Ideal Pulse Simulation](@ref) manual for details on the contents and format of the PET input file.

## 3. Environment settings

See [Field Simulation](@ref) manual for more details

## 4. Simulation settings

Simulation settings contain settings related to field simulation and pulse simulation. See [Field Simulation](@ref) and [Ideal Pulse Simulation](@ref) manuals for more details

## 5. DAQ settings

Current settings mimic the HADES DAQ setup. They contain the following fields:

- `"preamp"`: parameters modelling the preamplifier.
    - `"type"`: this field defines the rest of the parameters in the `"preamp"` field corresponding to given preamplifier type. Current available preamplifier type is `"generic"`, which defines a simple CR preamlifier. In the future other more complex models might be available. The following parameters in `"preamp"` describe a type `"generic"` preamplifier
    - `"t_decay_in_us"`: decay time constant in microseconds. Default 50 us.
    - `"t_rise_in_ns"`: rise time constant in nanoseconds. Default 100 ns.
    - `"gain_ADC_eV"`: preamplifier gain in units ADC/eV. Defined this way for convenience e.g. estimating gain based on the ADC value of a peak with known eV energy such as FEP etc. Default  0.138 ADC/eV.
    - `"offset_in_ADC"`: waveform baseline offset in units of ADC. Alternative input `"offset_in_keV"` to define offset in units of keV. Default 0.
    - `"noise_sigma_in_ADC"`: Gaus noise sigma in units of ADC. Alternative input `"noise_sigma_in_keV` to define Gaus noise sigma in units of keV. Default 0.

    WARNING: the implementation is janky. If you do not want to simulate noise, do NOT provide noise sigma in the preamplifier settings. Even if you provide noise settings for no noise (see below), if the sigma is non-zero, the noise WILL be added. WIP to fix this bad implementation.

- `"fadc"`: parameters modelling the FADC component of the DAQ chain.
    - `"type"`: this field defines the rest of the parameters in the `"fadc"` field corresponding to given FADC type.  Current available FADC type is `"generic"`, which defines a very simple dummy FADC component. In the future realistic FlashCam model might be available. The following parameters in `"fadc"` describe a type `"generic"` FADC
    - `"sampling_interval"`: sampling interval in nanoseconds. For example, in HADES data we see a time step of 16 ns. No default value
    - `"sampling_rate"`: alternative input to `"sampling_interval"`, sampling frequency of the FADC in MHz. No default value. One or the other parameters has to be provided.

    The `"generic"` type FADC is a very dummy model. It "samples" the waveform by simply picking every Nth point of the waveform based on the given sampling interval/frequency. For this you need to make sure that the time step of the simulated pulse is smaller than sampling interval. If the time step of the waveform does not fit into sampling interval in integer amount, the step will be rounded to int.
    The digitizing is simulated by converting the values to `Int16`. If the waveform ADC value surpasses maximum of `Int16`, the waveform will be saturated.

- `"trigger"`: parameteres modelling the trigger component of the DAQ chain.
    - `"type"`: this field defines the rest of the parameters in the `"trigger"` field corresponding to given Trigger type. Current available Trigger type is `"trapezoidal"`,  which defines a trapezoidal filter trigger. The following parameters in `"trigger"` describe a type `"trapezoidal"` Trigger
    - `"window_lengths"`: list of three trap filter window lengths in number of samples e.g. `[250,250,250]`. No default value.
    - `"threshold"`: threshold in keV. No default value.

    Pulses that don't pass the trigger threshold

- `"daq"`: parameters modelling data acquisition to disk.
    - `"type"`: this field defines the rest of the parameters in `"daq"` field corresponding to given DAQ type. Current available DAQ type is `"generic"`, which defines a simple dummy DAQ. The following parameters in `"daq"` describe a type `"generic"` DAQ
    - `"nsamples"`: number of samples in the waveform to save to disk. E.g. in HADES data there are 3748 samples in each waveform
    - `"baseline_length"`: number of samples of baseline to save before the trigger point.


## 6. Noise model

