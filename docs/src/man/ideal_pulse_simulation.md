# Ideal Pulse Simulation

As inputs to ideal pulse simulation you should provide

1. LEGEND detector metadata

2. PET input file

3. Environment settings

4. Simulation settings

5. Noise model

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

noise_model = Dict(
    "type" => "sim"
)

pss_table, pss_truth = LegendGeSim.simulate_pulses(detector_metadata_filename, pet_input_file, environment_settings, simulation_settings, noise_model; n_waveforms=100)
```

The optional variable `n_waveforms` defines how many waveforms will be simulted based on the given `pet` file. Not providing the variable will default `0` meaning, simulate all waveforms in the `pet` input file.

The output:
`pss_table` contains the following fields:
- `channel`: currently, ID of the contact (1 for p+ contact, 2 for n+ contact)
- `ievt`: event number
- `waveform`: waveform in format of `RDWaweform` from `RadiationDetectorSignals.jl`. Contains fields `time` and `signal`.

`pss_truth` contains a copy of the information in the `pet` input file for convenience (to compare)


## 1. LEGEND detector metadata

To obtain JSONs for all LEGEND detectors you can clone the repository `legend-detectors`. The JSONs can be found under `germanium/diodes`. You may also create your custom detector JSON file as long as it follows the LEGEND detector metadata format (see README in `diodes/`).

## 2. PET input file

The `pet` input file contains information on `p`osition, `e`nergy amd `t`ime of the depositions.

Current available formats:
- `g4simple` output in `hdf5` format. The `g4simple` simulations of HADES characterization data are available on MPIK under `/lfs/l1/legend/detector_char/enr/hades/simulations/legend-g4simple-simulation/simulations/`. Consult the `legend-g4simple-simulation` package to simulate your own custom `pet` input for given detector and source.

- `Geant4` output in `csv` format. See example in [legend-testdata](https://github.com/legend-exp/legend-testdata/) under `data/ldsim/single-invcoax-th228-geant4.csv`

Here we will use an example `pet` file `csv` format from `LegendTestData` corresponding to the Public Inverted Coax.

## 3. Environemnt settings

See [Field Simulation](@ref) manual for more details

## 4. Simulation settings

Simulation settings contain settings related to field simulation (see [Field Simulation](@ref) manual) + pulse simulation settings.

Options related to pulse simulation:

- `"time_step"`: time step in nanoseconds, default 1ns
- `"diffusion"`: setting for cloud charge simulation. Default `false` (point charges)
- `"self_repulsion"`: setting for cloud charge simulation. Default `false` (point charges)
- `"number_of_carriers"`: setting for cloud charge simulation. Default `1` (point charges)

## 5. Noise model

At this stage of the simulation, the noise model regards fano noise of the germanium crystal. The options are simply to simulate the noise or not.

If you would like to simulate fano noise, provide

```julia
noise_model = Dict( "type" => "sim" )
```

If you do not want to simulate fano noise, simply do not provide `noise_model` to to `simulate_pulses()` (it will default to `"type" => "none"`)

