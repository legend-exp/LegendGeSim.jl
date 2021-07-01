using SolidStateDetectors
using Plots
# plotlyjs(); # backend for interactive plots (Jupyter)

##
det_config_json = SSD_examples[:BEGe]
simulation = Simulation{Float32}(det_config_json)

## 
plot_det = plot(simulation.detector)
# png(plot_det, "plots_md/ssd_geometry.png")

##
using LegendGeSim

##
# import Pkg
# Pkg.add(url="https://github.com/legend-exp/LegendTestData.jl.git")
# Pkg.build("LegendTestData")

##
using LegendTestData
testdata_path = joinpath(LegendTestData.legend_test_data_path(), "data", "ldsim")


##
# path to my detector metadata
det_metadata = joinpath(testdata_path, "invcoax-metadata.json")
# convert LEGEND metadata to SSD configuration (only geometry!)
det_config_ssd = LegendGeSim.ssd_config(det_metadata)
# now we can plug it into SSD
simulation1 = Simulation(SolidStateDetector{Float32}(det_config_ssd))
# it's a bit boring cause we see the same detector as before
plot_leg = plot(simulation1.detector)
# png(plot_leg, "plots_md/metadata_geometry.png")


##
detector = LegendGeSim.simulate_detector(det_metadata, "configs/detector_study_ssd.json")

##
calculate_capacitance(detector)
get_active_volume(detector.point_types)

##
plot_ef = plot(detector.electric_field, φ = 0.0)
plot_electric_fieldlines!(detector, φ = 0.0)
# png(plot_ef, "plots_md/metadata_ef.png")

##
detector = LegendGeSim.simulate_detector(det_metadata, "configs/detector_study_fieldgen.json")