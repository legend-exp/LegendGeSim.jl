using LegendGeSim
using HDF5
using LegendDataTypes

##

mc_name = "raw-IC160A-Th228-uncollimated-top-run0002-source_holder-bi-hdf5-01-test"
mc_path = "data/"
processed_dir = "cache/"

## simulate waveforms with SSD, simulate noise from scratch
# sim_config_filename = "data/sim_config_SSD_NoiseSim.json"

## simulate waveforms with SSD, simulate noise and offset from data baselines
sim_config_filename = "data/sim_config_SSD_NoiseData.json"

sim_config = LegendGeSim.load_config(sim_config_filename)

filename(path) = splitext(basename(path))[1]
##

@info "----- g4simple -> mcstp"
g4_name = joinpath(mc_path, mc_name * ".hdf5")
mcstp_table = LegendGeSim.g4_to_mcstp(g4_name, sim_config.detector)

##
# preliminary step, you can save it to use later, or skip saving and use mcstp_table in memory for next steps
@info "Saving table"
det_name = filename(sim_config.detector)
out_filename = joinpath(processed_dir, "$(mc_name)_$(det_name)_mcstp.h5")
h5open(out_filename, "w") do f
    LegendDataTypes.writedata(f, "mcstp", mcstp_table)
end
@info "-> $out_filename"

##

# Read mcstp (if the previous step was skipped, otherwise skip this step)
@info "Reading mcstp"
det_name = filename(sim_config.detector)
mcstp_name = joinpath(processed_dir, "$(mc_name)_$(det_name)_mcstp.h5")
mcstp_table = HDF5.h5open(mcstp_name, "r") do input
    LegendDataTypes.readdata(input, "mcstp")
end

##
@info "----- mcstp -> mcpss"
mcpss_table, mcpss_mctruth = LegendGeSim.mcstp_to_mcpss(mcstp_table, sim_config_filename)

##
@info "Saving table"
det_name = filename(sim_config.detector)
out_filename = joinpath(processed_dir, "$(mc_name)_$(det_name)_$(sim_config.sim_method)_mcpss.h5")
h5open(out_filename, "w") do f
   LegendDataTypes.writedata(f, "mcpss/mcpss", mcpss_table)
   LegendDataTypes.writedata(f, "mcpss/mctruth", mcpss_mctruth)
end
@info "-> $out_filename"

##

# Read mcpss (if the previous step was skipped, otherwise skip this step)
det_name = filename(sim_config.detector)
mcpss_name = joinpath(processed_dir, "$(mc_name)_$(det_name)_$(sim_config.sim_method)_mcpss.h5")

mcpss_table = LegendGeSim.read_mcpss(mcpss_name)
mcpss_mctruth = LegendGeSim.read_mctruth(mcpss_name)

##
@info "----- mcpss -> mcraw"
mcraw_table = LegendGeSim.mcpss_to_mcraw(mcpss_table, mcpss_mctruth, sim_config_filename) 

##
@info "Saving table"
det_name = filename(sim_config.detector)
noise_name = occursin("Data", sim_config_filename) ? "NoiseData" : "NoiseSim"
out_filename = joinpath(processed_dir, "$(mc_name)_$(det_name)_$(sim_config.sim_method)_$(noise_name)_mcraw.h5")
HDF5.h5open(out_filename, "w") do f
   LegendDataTypes.writedata(f, "raw", mcraw_table)
   LegendDataTypes.writedata(f, "mctruth", mcpss_mctruth)
end
@info "-> $out_filename"
