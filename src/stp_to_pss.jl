function stp_to_pss(stp_table::Table, sim_config_file::AbstractString)
    # sim_config = load_config(sim_config_file)
    sim_config = propdict(sim_config_file)
    # basename for future cached files
    # TODO: rather than taking config name, ask user what name they want to use?
    # construct one based on detector, environment and PSS model?
    # (since noise model and elec are irrelevant)
    config_name = splitext(basename(sim_config_file))[1]
    stp_to_pss(stp_table, sim_config, config_name)
end


function stp_to_pss(stp_table::Table, sim_config::PropDict, config_name::AbstractString)
    det_meta = propdict(sim_config.detector_metadata)
    # det_meta = PropDicts.read(PropDict, sim_config.detector_metadata)
    env = Environment(sim_config)
    ps_simulator = PSSimulator(sim_config)
    noise_model = NoiseModel(sim_config)

    stp_to_pss(stp_table, det_meta, env, ps_simulator, noise_model, config_name)

end


"""
    stp_to_pss(stp_table, det_meta, env, ps_simulator, noise_model, config_name)

Table, PropDict, Env, PSSimulator, NoiseModel, AbstractString -> Table, Table    

Simulate waveforms based on stepping info given in <stp_table>, 
    LEGEND detector metadata <det_meta>, environment settings given in <env>,
    pulse shape simulation method <ps_simulator>, and <noise_model>

The output is a table with simulated pulses, and a table with simulation truth
    (may be abolished in the future since basically corresponds to stp table)
"""
function stp_to_pss(stp_table::Table, det_meta::PropDict, env::Environment, ps_simulator::PSSimulator, noise_model::NoiseModel,
    config_name::AbstractString)
    @info "---------------------- stp -> pss (pulse shape simulation)"

    # add fano noise, don't add if data noise is applied later
    # note: current understanding is that Siggen simulation does NOT include fano noise
    @info "//\\//\\//\\ Fano noise"
    stp_table_fano = fano_noise(stp_table, det_meta, env, noise_model)

    @info "_||_||_||_ Simulate detector" 
    detector = simulate_detector(det_meta, env, config_name, ps_simulator)

    @info "~.~.~.~.~ Simulate charge pulses"
    simulate_waveforms(stp_table_fano, detector)
end





