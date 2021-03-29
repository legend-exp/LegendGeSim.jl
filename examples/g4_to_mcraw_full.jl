using LegendGeSim
using HDF5
using LegendDataTypes

##

# full chain g4 -> mcstp -> mcpss -> mcraw
# only the last stage is saved

mc_name = "raw-IC160A-Th228-uncollimated-top-run0002-source_holder-bi-hdf5-01-test"
mc_path = "data/"
sim_config = "data/sim_datanoise.json"
processed_dir = "cache/"

##

mc_file = joinpath(mc_path, mc_name*".hdf5")
mcraw_table, mcpss_mctruth = LegendGeSim.g4_to_mcraw(mc_file, sim_config)

##
@info "Saving table"
out_filename = joinpath(processed_dir, mc_name*"_mcraw.h5")
if !ispath(dirname(out_filename)) mkpath(dirname(out_filename)) end
HDF5.h5open(out_filename, "w") do f
    LegendDataTypes.writedata(f, "raw", mcraw_table)
    LegendDataTypes.writedata(f, "mctruth", mcpss_mctruth)
end
@info("-> $out_filename")