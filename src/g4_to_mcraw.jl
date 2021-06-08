"""
    g4_to_mcraw(mc_file_path, sim_config_filename)

Simulate the full chain g4simple->mcstp->mcpss->mcraw

mc_file_path: full path to a g4simple file
sim_config_filename: simulation configuration json file
"""
function g4_to_mcraw(mc_file_path::AbstractString, sim_config_filename::AbstractString)
    sim_config = load_config(sim_config_filename)    

    g4_to_mcraw(mc_file_path, sim_config)
end


"""
    g4_to_mcraw(mc_file_path, sim_config)

Simulate the full chain g4simple->mcstp->mcpss->mcraw

mc_file_path: full path to a g4simple file
sim_config: PropDict object based on the simulation configuration json file
"""
function g4_to_mcraw(mc_file_path::AbstractString, sim_config::PropDict)

    @info "----- g4simple -> mcstp"
    mcstp_table = g4_to_mcstp(mc_file_path, sim_config.detector)

    @info "----- mcstp -> mcpss"
    ps_simulator = PSSimulator(sim_config)
    noise_model = NoiseModel(sim_config)
    mcpss_table, mcpss_mctruth = mcstp_to_mcpss(mcstp_table, sim_config.detector, ps_simulator, noise_model)

    @info "----- mcpss -> mcraw"
    daq = DAQmodel(sim_config)
    preamp = PreAmp(sim_config)
    mcraw_table = mcpss_to_mcraw(mcpss_table, mcpss_mctruth, daq, preamp, noise_model) 

    mcraw_table, mcpss_mctruth
end
