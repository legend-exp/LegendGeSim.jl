# Get detector metadata JSON from `legend-test-data`

using LegendTestData

ldsim_path = joinpath(legend_test_data_path(), "data", "ldsim")

detector_name = "invcoax-metadata"
detector_metadata_filename = joinpath(ldsim_path, detector_name*".json");

using LegendGeSim
using Plots

detector = LegendGeSim.LEGEND_SolidStateDetector(detector_metadata_filename)

plot(detector, aspect_ratio=1.0, camera=(20,20), size=(450,450), title=detector_name)

savefig("tutorial_det.pdf") # hide
