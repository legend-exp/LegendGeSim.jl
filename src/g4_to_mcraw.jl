
function g4_to_mcraw(det_name::AbstractString, det_path::AbstractString, mc_name::AbstractString, mc_path::AbstractString, processed_dir::AbstractString)

    @info "----- g4simple -> mcstp"
    mcstp_table = LegendGeSim.g4_to_mcstp(joinpath(mc_path, mc_name * ".hdf5"))

    @info "----- mcstp -> mcpss"
    mcpss_table, mcpss_mctruth = LegendGeSim.mcstp_to_mcpss(det_path, det_name, mcstp_table)

    @info "----- mcpss -> mcraw"
    mcraw_table = LegendGeSim.mcpss_to_mcraw(mcpss_table, mcpss_mctruth) # saves the final output

    mcraw_table, mcpss_mctruth
end
