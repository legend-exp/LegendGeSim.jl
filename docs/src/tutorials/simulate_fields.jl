using LegendGeSim
using Plots

using LegendTestData

germanium_testdata_path = joinpath(legend_test_data_path(), "data", "legend", "metadata", "hardware", "detectors", "germanium")

detector_name = "V99000A"
detector_metadata_filename = joinpath(germanium_testdata_path, "diodes", detector_name*".json");

environment_settings = Dict(
    "crystal_temperature_in_K" => 77,
    "medium" => "vacuum",
    "dl" => "vendor" # optional, default 0
);

simulation_settings = Dict(
    "method" => "SSD",
    "cached_name" => "test",
    "crystal_metadata_path" => joinpath(germanium_testdata_path, "crystals")
);

sim_ssd = LegendGeSim.simulate_fields(detector_metadata_filename, environment_settings, simulation_settings; overwrite=true)

LegendGeSim.capacitance_matrix(sim_ssd)

using LegendGeSim: SolidStateDetectors as SSD

SSD.is_depleted(sim_ssd.point_types)

SSD.estimate_depletion_voltage(sim_ssd)

plot(
    plot(sim_ssd.point_types, full_det = true),
    plot(sim_ssd.electric_potential, full_det = true),
    plot(sim_ssd.weighting_potentials[1], full_det = true),
    plot(sim_ssd.weighting_potentials[2], full_det = true),
    begin
        plot(sim_ssd.electric_field, full_det = true)
        SSD.plot_electric_fieldlines!(sim_ssd, full_det = true)
    end,
    size = (1000, 800), layout = (3, 2)
)

savefig("SSD_multiplot.pdf") # hide

simulation_settings_siggen = deepcopy(simulation_settings)
simulation_settings_siggen["method"] = "fieldgen"

sim_fieldgen = LegendGeSim.simulate_fields(detector_metadata_filename, environment_settings, simulation_settings_siggen; overwrite=true)

LegendGeSim.capacitance_matrix(sim_fieldgen)

simulation_settings_siggen1 = Dict(
    "method" => "fieldgen",
    "cached_name" => "test",
    "impurity_profile" => "cache/V99000.dat",
    "offset_in_mm" => 0
);

sim_fieldgen = LegendGeSim.simulate_fields(detector_metadata_filename, environment_settings, simulation_settings_siggen1; overwrite=true)

e_pot = SSD.ElectricPotential(sim_fieldgen);
w_pot = SSD.WeightingPotential(sim_fieldgen);

plot(
    plot(sim_ssd.electric_potential, full_det = true, title = "Electric potential via SSD"),
    plot(e_pot, full_det = true, title = "Electric potential via Fieldgen"),
    size = (1000, 500), layout = (1, 2)
)

savefig("SSD_vs_fieldgen_Epot.pdf") # hide

plot(
    plot(sim_ssd.weighting_potentials[1], full_det = true, title = "Weighting potential via SSD"),
    plot(w_pot, full_det = true, title = "Weighting potential via fieldgen"),
    size = (1000, 500), layout = (1, 2)
)

savefig("SSD_vs_fieldgen_Wpot.pdf") # hide
