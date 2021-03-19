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
function g4_to_mcraw(det_name::AbstractString, det_path::AbstractString, mc_name::AbstractString, mc_path::AbstractString, sim_config::AbstractString)

    @info "----- g4simple -> mcstp"
    mcstp_table = LegendGeSim.g4_to_mcstp(joinpath(mc_path, mc_name * ".hdf5"))

    daq = construct_GenericDAQ(sim_config)
    preamp = construct_PreAmp(sim_config)
    noise_model = NoiseModel(sim_config)

    @info "----- mcstp -> mcpss"
    mcpss_table, mcpss_mctruth = LegendGeSim.mcstp_to_mcpss(det_path, det_name, mcstp_table, noise_model)

    @info "----- mcpss -> mcraw"
    mcraw_table = LegendGeSim.mcpss_to_mcraw(mcpss_table, mcpss_mctruth, daq, preamp, noise_model) 
    # mcraw_table = LegendGeSim.mcpss_to_mcraw(mcpss_table, mcpss_mctruth, sim_config) 

    mcraw_table, mcpss_mctruth
end
