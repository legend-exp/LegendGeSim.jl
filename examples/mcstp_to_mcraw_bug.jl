using LegendGeSim
using LegendDataTypes
using HDF5

mc_name = "raw-IC160A-Th228-uncollimated-top-run0002-source_holder-bi-hdf5-02"
det_name = "V05266A"
det_path = "data/"
processed_dir = "output/"
##

@info "Reading mcstp"
mcstp_name = joinpath("data", mc_name*"_mcstp.h5")
mcstp_table = HDF5.h5open(mcstp_name, "r") do input
    LegendDataTypes.readdata(input, "mcstp")
end

@info "----- mcstp -> mcpss"
mcpss_table, mcpss_mctruth = LegendGeSim.mcstp_to_mcpss(det_path, det_name, mcstp_table)

@info "----- mcpss -> mcraw"
# the bug is here
mcraw_table = LegendGeSim.mcpss_to_mcraw(mcpss_table, mcpss_mctruth) 

@info "Saving table..."
out_filename = joinpath(processed_dir, mc_name*"_mcraw.h5")
@info("-> $out_filename")
HDF5.h5open(out_filename, "w") do f
    LegendDataTypes.writedata(f, "raw", mcraw_table)
    LegendDataTypes.writedata(f, "mctruth", mcpss_mctruth)
end
@info("Done")
