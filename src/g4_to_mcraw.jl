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
    mcstp_table = LegendGeSim.g4_to_mcstp(mc_file_path)

    @info "----- mcstp -> mcpss"
    noise_model = LegendGeSim.NoiseModel(sim_config)
    simulation = LegendGeSim.detector_simulation(sim_config.detector_path, sim_config.detector)
    mcpss_table, mcpss_mctruth = LegendGeSim.mcstp_to_mcpss(simulation, mcstp_table, noise_model)

    @info "----- mcpss -> mcraw"
    daq = GenericDAQ(sim_config)
    preamp = PreAmp(sim_config)
    mcraw_table = LegendGeSim.mcpss_to_mcraw(mcpss_table, mcpss_mctruth, daq, preamp, noise_model) 

    mcraw_table, mcpss_mctruth
end
