# This file is a part of LegendGeSim.jl, licensed under the MIT License (MIT).

using SolidStateDetectors
using Plots
using LegendGeSim
using LegendTestData
using Unitful

testdata_path = joinpath(LegendTestData.legend_test_data_path(), "data", "ldsim")

detector_metadata_file = joinpath(testdata_path, "invcoax-metadata.json") # path to my detector metadata

################################
### using SolidStateDetectors.jl

simulation_config_file = joinpath(dirname(dirname(pathof(LegendGeSim))), "examples/configs/SSD_NoiseSim.json")
fields_sim_file = "det_fields_ssd" # file prefix for automated storing of simulated fields
geant4_output_hit_file = joinpath(testdata_path, "single-invcoax-th228-geant4.csv") # Geant 4 output hits

simulation_settings = LegendGeSim.SSDSimulator(LegendGeSim.propdict(simulation_config_file))
environment_settings = LegendGeSim.Environment(LegendGeSim.propdict(simulation_config_file))

sim = LegendGeSim.simulate_fields(detector_metadata_file, environment_settings, simulation_settings, fields_sim_file; overwrite = true)

plot(
    plot(sim.detector),
    plot(sim.point_types, full_det = true),
    plot(sim.electric_potential, full_det = true),
    begin
        plot(sim.electric_field, full_det = true, φ = 0.0)
        plot_electric_fieldlines!(sim, full_det = true, φ = 0.0)
    end,
    plot(sim.weighting_potentials[1], full_det = true),
    plot(sim.weighting_potentials[2], full_det = true),
    size = (1000, 1400),
    layout = (3, 2),
    full_det = true
)

C = calculate_capacitance_matrix(sim)
active_volume = get_active_volume(sim.point_types)

begin
    geant4_output_evt_table = LegendGeSim.read_pet(geant4_output_hit_file)
    # preprocessing of geant4 output: `pet` -> `stp` format
    stp_table = LegendGeSim.pet_to_stp(geant4_output_evt_table, sim)

    pss_table, pss_truth_table = LegendGeSim.stp_to_pss(stp_table, sim, simulation_settings)

    elec_chain = LegendGeSim.ElecChain(LegendGeSim.propdict(simulation_config_file))
    trigger = LegendGeSim.Trigger(LegendGeSim.propdict(simulation_config_file))
    daq = LegendGeSim.DAQ(LegendGeSim.propdict(simulation_config_file))
    noise_model = LegendGeSim.NoiseModel(LegendGeSim.propdict(simulation_config_file))

    raw_table = LegendGeSim.pss_to_raw(pss_table, pss_truth_table, simulation_settings, elec_chain, trigger, daq, noise_model)
end
plot(point_contact_waveforms.waveform, legend = false)
point_contact_waveforms = filter(row -> row.channel == 1 && row.energy > 5u"keV", raw_table)

raw_table_ssd = raw_table # For later comparison to SigGen




#########################
### using Field- & SigGen

simulation_config_file = joinpath(dirname(dirname(pathof(LegendGeSim))), "examples/configs/siggen_NoiseSim.json")
fields_sim_file = "det_fields" # file prefix for automated storing of simulated fields
geant4_output_hit_file = joinpath(testdata_path, "single-invcoax-th228-geant4.csv") # Geant 4 output hits

simulation_settings = LegendGeSim.SiggenSimulator(LegendGeSim.propdict(simulation_config_file))
environment_settings = LegendGeSim.Environment(LegendGeSim.propdict(simulation_config_file))

cd("LegendGeSim.jl/examples") # This has to be solved...
sim = LegendGeSim.simulate_fields(detector_metadata_file, environment_settings, simulation_settings, fields_sim_file; overwrite = true)

begin
    geant4_output_evt_table = LegendGeSim.read_pet(geant4_output_hit_file)
    #     # preprocessing of geant4 output: `pet` -> `stp` format
    stp_table = LegendGeSim.pet_to_stp(geant4_output_evt_table, sim)

    pss_table, pss_truth_table = LegendGeSim.stp_to_pss(stp_table, sim, simulation_settings)

    elec_chain = LegendGeSim.ElecChain(LegendGeSim.propdict(simulation_config_file))
    trigger = LegendGeSim.Trigger(LegendGeSim.propdict(simulation_config_file))
    daq = LegendGeSim.DAQ(LegendGeSim.propdict(simulation_config_file))
    noise_model = LegendGeSim.NoiseModel(LegendGeSim.propdict(simulation_config_file))

    raw_table = LegendGeSim.pss_to_raw(pss_table, pss_truth_table, simulation_settings, elec_chain, trigger, daq, noise_model)
end
point_contact_waveforms = filter(row -> row.channel == 1 && row.energy > 5u"keV", raw_table)
plot(point_contact_waveforms.waveform, legend = false)

raw_table_siggen = raw_table

################

begin
    i = 3
    wvssd = raw_table_ssd.waveform[i]
    wvsiggen = raw_table_siggen.waveform[i]
    plot(wvssd, label = "SSD")
    plot!(wvsiggen, label = "SigGen")
end
