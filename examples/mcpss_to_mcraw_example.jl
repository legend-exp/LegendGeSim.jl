using LegendGeSim
using LegendDataTypes
using HDF5
##

## Use this to debug mcpss->mcraw and tune DAQ and preamp parameters
# Reads an already existing mcpss file (to save time)

##

mc_name = "raw-IC160A-Th228-uncollimated-top-run0002-source_holder-bi-hdf5-01-test"
# mc_name = "raw-IC160A-Th228-uncollimated-top-run0002-source_holder-bi-hdf5-02"
mc_path = "cache/"
det_name = "V05266A"
det_path = "data/"
processed_dir = "cache/"
# daq_config = 'data/daq_config.json'
# preamp_config = 'data/preamp_config.json'
elec_config = "data/elec.json"

##
mcpss_name = joinpath(mc_path, mc_name*"_mcpss.h5")

# mcpss_table = read_mcpss(mcpss_name)
# mcpss_mctruth = read_mctruth(mcpss_name)
mcpss_table = LegendGeSim.read_mcpss(mcpss_name)
mcpss_mctruth = LegendGeSim.read_mctruth(mcpss_name)

##

@info "----- mcpss -> mcraw"
mcraw_table = LegendGeSim.mcpss_to_mcraw(mcpss_table, mcpss_mctruth, elec_config) 

##

@info "Saving table..."
out_filename = joinpath(processed_dir, mc_name*"_mcraw.h5")
@info("-> $out_filename")
HDF5.h5open(out_filename, "w") do f
    LegendDataTypes.writedata(f, "raw", mcraw_table)
    LegendDataTypes.writedata(f, "mctruth", mcpss_mctruth)
end
@info("Done")
