## Example usage of LegendGeSim package
## Using the function g4_to_mcstp, process a g4simple output into mcstp format
# (input to SSD waveform simulation)
# In this code, an example g4simple hdf5 is read,
# and the resulting mcstp table is saved to cache/

using LegendGeSim

using HDF5
using LegendDataTypes

##
# g4_dir = "/lfs/l1/legend/detector_char/enr/hades/simulations/legend-g4simple-simulation/IC-legend/IC160A/Th228/uncollimated/top_source_holder/hdf5/"
g4_dir = "data"
processed_dir = "cache"
base_filename = "raw-IC160A-Th228-uncollimated-top-run0002-source_holder-bi-hdf5-01-test"
raw_extension = ".hdf5"
processed_extension = ".h5"
##
@info "Processing g4 events"
full_name = joinpath(g4_dir, base_filename) * raw_extension
mcstp_table = LegendGeSim.g4_to_mcstp(full_name)

##

@info "Saving..."
out_filename = joinpath(processed_dir, base_filename) * processed_extension
if !ispath(dirname(out_filename)) mkpath(dirname(out_filename)) end

HDF5.h5open(out_filename, "w") do f
    LegendDataTypes.writedata(f, "mctruth", mcstp_table)
end
println("Processed file save to: $out_filename")
