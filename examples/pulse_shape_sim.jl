using LegendGeSim 

##
# import Pkg
# Pkg.add(url="https://github.com/legend-exp/LegendTestData.jl.git")
# Pkg.build("LegendTestData")
using LegendTestData
testdata_path = joinpath(LegendTestData.legend_test_data_path(), "data", "ldsim")

##

# PET (position-energy-time) information of the simulated energy depositions
pet_file = joinpath(testdata_path, "single-invcoax-th228-geant4.csv")
# detector
det_metadata = joinpath(testdata_path, "invcoax-metadata.json")
# simulation settings
# sim_config_name = "configs/SSD_NoiseSim.json"
sim_config_name = "configs/siggen_NoiseSim.json"


##
raw_table = LegendGeSim.simulate_raw(pet_file, det_metadata, sim_config_name)

##
using HDF5
using LegendHDF5IO

##
mkdir("output")
HDF5.h5open("output/my_name_raw.h5", "w") do f
    LegendHDF5IO.writedata(f, "raw", raw_table)
end

##
using Plots 

##
plot_wf = plot(raw_table.waveform[1:5])
# png(plot_wf, "raw.png")

##
stp_table = LegendGeSim.pet_to_stp(pet_file, det_metadata)



##
mkdir("cache")
HDF5.h5open("cache/my_name_stp.h5", "w") do f
    LegendHDF5IO.writedata(f, "stp", stp_table)
end

##
pss_table, pss_truth = LegendGeSim.stp_to_pss("cache/my_name_stp.h5", det_metadata, sim_config_name)

##
pss_table, pss_truth = LegendGeSim.stp_to_pss(stp_table, det_metadata, sim_config_name)

##
using Plots

##
plot_pss = plot(pss_table.waveform[1:5])
# png(plot_pss, "plots_md/pss.png")

##
h5open("cache/my_name_pss.h5", "w") do f
    LegendHDF5IO.writedata(f, "pss/pss", pss_table)
    LegendHDF5IO.writedata(f, "pss/truth", pss_truth)
end

##
raw_table = LegendGeSim.pss_to_raw("cache/my_name_pss.h5", det_metadata, sim_config_name)

##
raw_table = LegendGeSim.pss_to_raw(pss_table, pss_truth, det_metadata, sim_config_name)

##
pss_table, pss_truth = LegendGeSim.pet_to_pss(pet_filename, det_metadata, sim_config_name)
