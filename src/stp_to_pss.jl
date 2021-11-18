"""
    stp_to_pss(stp_table, det_meta, env, ps_simulator, noise_model, config_name)

Table, PropDict, Env, PSSimulator, NoiseModel, AbstractString -> Table, Table    

Simulate waveforms based on stepping info given in <stp_table>, 
    LEGEND detector metadata <det_meta>, environment settings given in <env>,
    pulse shape simulation method <ps_simulator>, and <noise_model>

The output is a table with simulated pulses, and a table with simulation truth
    (may be abolished in the future since basically corresponds to stp table)
"""
function stp_to_pss(stp_table::Table, det_meta::PropDict, env::Environment, ps_simulator::SSDSimulator, noise_model::NoiseModel,
    cached_name::AbstractString)

    @info "_||_||_||_ Simulate detector" 
    sim = simulate_fields(det_meta, env, ps_simulator, cached_name)

    @info "~.~.~.~.~ Simulate charge pulses"
    simulate_waveforms(stp_table, sim)
end


function stp_to_pss(stp_filepath::AbstractString, det_meta_fullpath::AbstractString, sim_config_filename::AbstractString)
    stp_table = h5open(stp_filepath, "r") do input
        LegendHDF5IO.readdata(input, "stp")
    end
    
    stp_to_pss(stp_table, det_meta_fullpath, sim_config_filename)
end

function stp_to_pss(stp_table::Table, sim::Union{Simulation, SigGenSetup}, sim_settings::PSSimulator)
    # Right now, only use `sim_settings` to dispatch for SSD 
    # Later, waveform simulation settings could be set in the configuration file and be parsed to `sim_settings`
    @info "~.~.~.~.~ Simulate charge pulses"
    simulate_waveforms(stp_table, sim)
end



function stp_to_pss(stp_table::Table, det_meta_fullpath::AbstractString, sim_config_filename::AbstractString)
    # construct simulation config based on given inputs and settings 
    sim_config = load_config(stp_table, det_meta_fullpath, sim_config_filename)

    stp_to_pss(sim_config)
end    


function stp_tp_pss(sim_config_filename::AbstractString)
    # given config file already contains inputs as well
    sim_config = propdict(sim_config_filename)
    # ok this looks kinda stupid i'm in a multiple dispatch loop what do i do
    stp_to_pss(sim_config.input_file, sim_config.detector_metadata, sim_config_filename)
end


function stp_to_pss(sim_config::PropDict)
    @info "---------------------- stp -> pss (pulse shape simulation)"

    det_meta = propdict(sim_config.detector_metadata)
    env = Environment(sim_config)
    ps_simulator = PSSimulator(sim_config)
    noise_model = NoiseModel(sim_config)

    stp_to_pss(sim_config.input_file, det_meta, env, ps_simulator, noise_model, sim_config.pss.cached_name)
end
