using LegendGeSim

using HDF5
using LegendDataTypes

##
det_name = "V05266A"
det_path = "data/"
mc_name = "raw-IC160A-Th228-uncollimated-top-run0002-source_holder-bi-hdf5-01-test"
mc_path = "data/"
processed_dir = "cache"

##
g4_name = joinpath(mc_path, mc_name * ".hdf5")
mcstp_table = LegendGeSim.g4_to_mcstp(g4_name)
##
mcpss_table, mcpss_mctruth = LegendGeSim.mcstp_to_mcpss(det_path, det_name, mcstp_table)

## save results
println("...saving table")
out_filename = joinpath(processed_dir, "$(mc_name)_mcpss.h5")
println("-> $out_filename")
HDF5.h5open(out_filename, "w") do f
    LegendDataTypes.writedata(f, "mcpss/mcpss", mcpss_table)
    LegendDataTypes.writedata(f, "mcpss/mctruth", mcpss_mctruth)
end
println("Done")
