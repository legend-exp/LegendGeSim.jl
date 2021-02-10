using LegendGeSim

using HDF5
using LegendDataTypes

##

det_name = "V05266A"
det_path = "data/"
mc_name = "raw-IC160A-Th228-uncollimated-top-run0002-source_holder-bi-hdf5-01-test"
mc_path = "data/"

# processed_dir = "cache"
processed_dir = "/home/sagitta/_legend/pss/LegendGeSim.jl/examples/cache/"

##

mcraw_table, mcpss_mctruth = LegendGeSim.g4_to_mcraw(det_name, det_path, mc_name, mc_path)

##
@info "Saving table..."
out_filename = joinpath(processed_dir, mc_name*"_mcraw.h5")
@info("-> $out_filename")
HDF5.h5open(out_filename, "w") do f
    LegendDataTypes.writedata(f, "raw", mcraw_table)
    LegendDataTypes.writedata(f, "mctruth", mcpss_mctruth)
end
@info("Done")
