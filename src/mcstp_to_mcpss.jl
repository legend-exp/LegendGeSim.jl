# Random.seed!(123) # only for testing


"""
    mstp_to_mcpss(mc_events, sim_config_file)

Simulate waveforms based on given MC events for a given configuration.
Returns the resulting simulated waveforms and MC truth

mc_events: table in the mcstp format (output of g4_to_mcstp)
sim_config_file: simulation configuration json file
"""
function mcstp_to_mcpss(mc_events::Table, sim_config_file::AbstractString)
    sim_config = LegendGeSim.load_config(sim_config_file)
    mcstp_to_mcpss(mc_events, sim_config)

end


function mcstp_to_mcpss(mc_events::Table, sim_config::PropDict)
    ps_simulator = PSSimulator(sim_config)
    @info("Reading geometry for $(sim_config.detector)")
    noise_model = NoiseModel(sim_config)

    mcstp_to_mcpss(mc_events, sim_config.detector, ps_simulator, noise_model)

end


"""
    mstp_to_mcpss(mc_events, sim_config)

Simulate waveforms based on given MC events for a given configuration.
Returns the resulting simulated waveforms and MC truth

mc_events: table in the mcstp format (output of g4_to_mcstp)
sim_config: PropDict object with simulation configuration
"""
function mcstp_to_mcpss(mc_events::Table, det_json::AbstractString, ps_simulator::PSSimulator, noise_model::NoiseModel)
    # add fano noise, don't add if data noise is applied later
    # should we do this before siggen as well? or does siggen do it by itself?
    # noise_model = NoiseModel(sim_config)
    mc_events = fano_noise(mc_events, det_json, noise_model)

    # mcpss_table, mcpss_mctruth = simulate_wf(mc_events, detector_config, ps_simulator)
    simulate_wf(mc_events, det_json, ps_simulator)

end





