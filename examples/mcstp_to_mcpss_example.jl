using LegendGeSim
using LegendDataTypes
using HDF5

##

## Use to save time to do only one stage mcstp->mcpss
# Reads an already existing mcstp file

mc_name = "raw-IC160A-Th228-uncollimated-top-run0002-source_holder-bi-hdf5-01-test"
mc_path = "cache/"
det_name = "IC160A"
det_path = "data/"
processed_dir = "cache/"
sim_config_file = "data/sim_datanoise.json"

##

@info "Reading mcstp"
mcstp_name = joinpath(mc_path, mc_name*"_mcstp.h5")
mcstp_table = HDF5.h5open(mcstp_name, "r") do input
    LegendDataTypes.readdata(input, "mcstp")
end

##

@info "----- mcstp -> mcpss"
sim_config = LegendGeSim.load_config(sim_config_file)
mcpss_table, mcpss_mctruth = LegendGeSim.mcstp_to_mcpss(det_path, det_name, mcstp_table, sim_config)

##

## save results
println("...saving table")
out_filename = joinpath(processed_dir, "$(mc_name)_mcpss.h5")
println("-> $out_filename")
HDF5.h5open(out_filename, "w") do f
    LegendDataTypes.writedata(f, "mcpss/mcpss", mcpss_table)
    LegendDataTypes.writedata(f, "mcpss/mctruth", mcpss_mctruth)
end
println("Done")
