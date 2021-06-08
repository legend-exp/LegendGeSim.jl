function mcstp_to_mcpss(mc_events::Table, sim_config_file::AbstractString)
    sim_config = load_config(sim_config_file)
    mcstp_to_mcpss(mc_events, sim_config)
end


"""
    mcstp_to_mcpss(mc_events, sim_config)

Simulate waveforms based on given MC events for a given configuration.
Returns the resulting simulated waveforms and MC truth

mc_events: table in the mcstp format (output of g4_to_mcstp)
sim_config: PropDict with the given JSON configuration
"""
function mcstp_to_mcpss(mc_events::Table, sim_config::PropDict)
    ps_simulator = PSSimulator(sim_config)
    noise_model = NoiseModel(sim_config)
    env = Environment(sim_config)

    mcstp_to_mcpss(mc_events, sim_config.detector, env, ps_simulator, noise_model)

end


"""
    mcstp_to_mcpss(mc_events, det_json, env, ps_simulator, noise_model)

Simulate waveforms based on given MC events for given
    - detector geometry
    - environment settings
    - pulse shape simulation method
    - noise model

Returns
    - table with resulting simulated waveforms in the mcpss format
    - table with MC truth

mc_events: table in the mcstp format (output of g4_to_mcstp)
det_json: path to detector LEGEND metadata JSON 
env: Environment object with parameters user defined in simulation config
"""
function mcstp_to_mcpss(mc_events::Table, det_json::AbstractString, env::Environment, ps_simulator::PSSimulator, noise_model::NoiseModel)
    # add fano noise, don't add if data noise is applied later
    # should we do this before siggen as well? or does siggen do it by itself?
    # noise_model = NoiseModel(sim_config)
    mc_events = fano_noise(mc_events, det_json, env, noise_model)

    simulate_wf(mc_events, det_json, env, ps_simulator)
end





