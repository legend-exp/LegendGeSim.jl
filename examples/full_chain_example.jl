using LegendGeSim
using LegendTestData
using LegendGeSim: SolidStateDetectors as SSD
using Plots
using Unitful

# We need 3 input files:
# 1) Detector meta data configuration file
# 2) Energy depositions hits (e.g. output from G4Simple)
# 3) LegendGeSim-configuration file

testdata_path = joinpath(LegendTestData.legend_test_data_path(), "data", "ldsim")
LegendGeSimExampleConfigFolder = joinpath(dirname(dirname(pathof(LegendGeSim))), "examples/configs")

detector_metadata_filename = "invcoax-metadata.json"
geant4_output_hit_filename = "single-invcoax-th228-geant4.csv"
cp(joinpath(testdata_path, detector_metadata_filename), detector_metadata_filename, force = true)
cp(joinpath(testdata_path, geant4_output_hit_filename), geant4_output_hit_filename, force = true)

sim_settings_ssd_filename = "SSD_NoiseSim.json"
sim_settings_siggen_filename = "siggen_NoiseSim.json"
cp(joinpath(LegendGeSimExampleConfigFolder, sim_settings_ssd_filename), sim_settings_ssd_filename, force = true)
cp(joinpath(LegendGeSimExampleConfigFolder, sim_settings_siggen_filename), sim_settings_siggen_filename, force = true)

## Load the config
config = LegendGeSim.load_config(detector_metadata_filename, sim_settings_ssd_filename);
# config = LegendGeSim.load_config(detector_metadata_filename, sim_settings_siggen_filename);

## Simulate the fields of the detector
sim = LegendGeSim.simulate_fields(config, overwrite = true);

plot(
    plot(sim.detector),
    plot(sim.point_types, full_det = true),
    plot(sim.electric_potential, full_det = true),
    begin
        plot(sim.electric_field, full_det = true, φ = 0.0)
        SSD.plot_electric_fieldlines!(sim, full_det = true, φ = 0.0)
    end,
    plot(sim.weighting_potentials[1], full_det = true),
    plot(sim.weighting_potentials[2], full_det = true),
    size = (1000, 1400),
    layout = (3, 2),
    full_det = true
)

C = LegendGeSim.capacitance_matrix(sim)


## Load the energy depositions hits 
geant4_output_evt_table = LegendGeSim.read_pet(geant4_output_hit_filename); # could also in the HDF5 Format -> LegendTextIO.jl & LegendHDF5IO.jl

### Preprocessing of the events. 
# In case of SSD: Filters depositions not located inside the detector
stp_table = LegendGeSim.pet_to_stp(geant4_output_evt_table, sim);

### Simulation of charge drifts through the detector and generate signals (waveforms) either via SSD or via Siggen.
# These are "perfect": No noise.
pss_table, pss_truth_table = LegendGeSim.stp_to_pss(stp_table, sim, config);

### Generate Data-Like Raw waveforms
# Add Noise, Electronic Response (from PreAmp), simulate DAQ
raw_table = LegendGeSim.pss_to_raw(pss_table, pss_truth_table, config)

plot(raw_table.waveform[1])
