"""
    g4_to_mcraw(det_name, det_path, mc_name, mc_path)

Simulate the full chain g4->mcstp->mcpss->mcraw

det_name: detector name (e.g. "V05266A")
det_path: path to detector json file
mc_name: name of the g4simple hdf5 file
mc_path: path to the folder containing the g4simple hdf5 file
elec_config: path to elec config json file 

Returns two tables: mcraw and mctruth information
"""

function g4_to_mcraw(mc_file_path::AbstractString, sim_config_filename::AbstractString)
    sim_config = load_config(sim_config_filename)    

    g4_to_mcraw(mc_file_path, sim_config)
end


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
