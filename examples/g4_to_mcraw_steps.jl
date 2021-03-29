using LegendGeSim
using HDF5
using LegendDataTypes

##

mc_name = "raw-IC160A-Th228-uncollimated-top-run0002-source_holder-bi-hdf5-01-test"
mc_path = "data/"
processed_dir = "cache/"
sim_config_filename = "data/sim_datanoise.json"

##

@info "----- g4simple -> mcstp"
mcstp_table = LegendGeSim.g4_to_mcstp(joinpath(mc_path, mc_name * ".hdf5"))

##
# preliminary step, you can save it to use later, or skip saving and use mcstp_table in memory for next steps
@info "Saving table"
out_filename = joinpath(processed_dir, "$(mc_name)_mcstp.h5")
h5open(out_filename, "w") do f
    LegendDataTypes.writedata(f, "mcstp", mcstp_table)
end
@info "-> $out_filename"

##

# Read mcstp (if the previous step was skipped, otherwise skip this step)
@info "Reading mcstp"
mcstp_name = joinpath(processed_dir, mc_name*"_mcstp.h5")
mcstp_table = HDF5.h5open(mcstp_name, "r") do input
    LegendDataTypes.readdata(input, "mcstp")
end

# Choose version 1 or version 2 (equivalent)
##
# version 1
@info "----- mcstp -> mcpss"
mcpss_table, mcpss_mctruth = LegendGeSim.mcstp_to_mcpss(mcstp_table, sim_config_filename)

##
# version 2
@info "----- mcstp -> mcpss"
sim_config = LegendGeSim.load_config(sim_config_filename)
mcpss_table, mcpss_mctruth = LegendGeSim.mcstp_to_mcpss(mcstp_table, sim_config)

##
# version 3
sim_config = LegendGeSim.load_config(sim_config_filename)
simulation = LegendGeSim.detector_simulation(sim_config.detector_path, sim_config.detector)
ps_simulator = LegendGeSim.PSSimulator(sim_config)
noise_model = LegendGeSim.NoiseModel(sim_config)
mcpss_table, mcpss_mctruth = LegendGeSim.mcstp_to_mcpss(mcstp_table, simulation, ps_simulator, noise_model)

##

@info "Saving table"
out_filename = joinpath(processed_dir, "$(mc_name)_mcpss.h5")
h5open(out_filename, "w") do f
   LegendDataTypes.writedata(f, "mcpss/mcpss", mcpss_table)
   LegendDataTypes.writedata(f, "mcpss/mctruth", mcpss_mctruth)
end
@info "-> $out_filename"

##

# Read mcpss (if the previous step was skipped, otherwise skip this step)
mcpss_name = joinpath(processed_dir, mc_name*"_mcpss.h5")

mcpss_table = LegendGeSim.read_mcpss(mcpss_name)
mcpss_mctruth = LegendGeSim.read_mctruth(mcpss_name)

# Choose version 1, 2 or 3 (equivalent)
##
# version 1
# mcpss_to_mcraw will load the simulation config
@info "----- mcpss -> mcraw"
mcraw_table = LegendGeSim.mcpss_to_mcraw(mcpss_table, mcpss_mctruth, sim_config_filename) 

##
# version 2
# we already loaded the simulation config before, can pass it to mcpss_to_mcraw
@info "----- mcpss -> mcraw"
mcraw_table = LegendGeSim.mcpss_to_mcraw(mcpss_table, mcpss_mctruth, sim_config) 

##
# version 3
# we already loaded the simulation config and constructed noise_model before, pass both to mcpss_to_mcraw
@info "----- mcpss -> mcraw"
daq = LegendGeSim.GenericDAQ(sim_config)
preamp = LegendGeSim.PreAmp(sim_config)
mcraw_table = LegendGeSim.mcpss_to_mcraw(mcpss_table, mcpss_mctruth, daq, preamp, noise_model) 


##
@info "Saving table"
out_filename = joinpath(processed_dir, mc_name*"_mcraw.h5")
HDF5.h5open(out_filename, "w") do f
   LegendDataTypes.writedata(f, "raw", mcraw_table)
   LegendDataTypes.writedata(f, "mctruth", mcpss_mctruth)
end
@info "-> $out_filename"
