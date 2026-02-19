# # Simulate Fields

#md # You can also download this tutorial as a
#md # [Jupyter notebook](simulate_fields.ipynb) and a plain
#md # [Julia source file](simulate_fields.jl).

# To obtain JSONs for all LEGEND detectors you can clone the repository `legend-detectors`. The JSONs can be found under `germanium/diodes`. You may also create your custom detector JSON file as long as it follows the LEGEND detector metadata format (see README in `diodes/`).
#
# Here we will use public example JSONs from `LegendTestData`

using LegendGeSim
using Plots

# ## Get detector metadata JSON from `legend-test-data`

using LegendTestData
#
germanium_testdata_path = joinpath(legend_test_data_path(), "data", "legend", "metadata", "hardware", "detectors", "germanium")

detector_name = "V99000A"
detector_metadata_filename = joinpath(germanium_testdata_path, "diodes", detector_name*".json");

# Alternatively, enter your own path to a real LEGEND detector JSON (see manual on Field Simulation for more details)
#
# ```julia
# detector_metadata_filename = "path/to/V04545A.json"
# ```

# ## Settings

# See manual on [Field Simulation](@ref) for a detailed explanation of environment and simulation settings

environment_settings = Dict(
    "crystal_temperature_in_K" => 77, 
    "medium" => "vacuum", 
    "dl" => "vendor" # optional, default 0
);

#

simulation_settings = Dict(
    "method" => "SSD",
    "cached_name" => "", # a non-empty string will cache the simulation results
    "crystal_metadata_path" => joinpath(germanium_testdata_path, "crystals")
);

# ## Simulation with `SolidStateDetectors`

sim_ssd = LegendGeSim.simulate_fields(detector_metadata_filename, environment_settings, simulation_settings; overwrite=true)

# The output of `simulate_fields()` in case of SSD will be a `SolidStateDetector` object.
#
# Below are some examples of what you can do with it

LegendGeSim.capacitance_matrix(sim_ssd)

# Using built-in `SolidStateDetectors.jl` functionality

using LegendGeSim: SolidStateDetectors as SSD

#

SSD.is_depleted(sim_ssd.point_types)

#

SSD.estimate_depletion_voltage(sim_ssd)

#

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

#jl savefig("SSD_multiplot.pdf") # hide
#md savefig("SSD_multiplot.svg"); nothing # hide
#md # ![SSD_multiplot](SSD_multiplot.svg)


# Consult `SolidStateDetectors.jl` tutorials for more operations on the `SolidStateDetector` object

# ## Simulation with `Fieldgen`

# The same simulation settings can be used for siggen, changing only the `"method"` settings to `"siggen"` (or `"fieldgen"`) and loading the `MJDSigGen` package.

using MJDSigGen
simulation_settings_siggen = deepcopy(simulation_settings)
simulation_settings_siggen["method"] = "fieldgen"

# The same function `simulate_fields` is used with identical inputs as in case of `SSD`.

sim_fieldgen = LegendGeSim.simulate_fields(detector_metadata_filename, environment_settings, simulation_settings_siggen; overwrite=true)

# The output of `simulate_fields()` in case of Fieldgen will be a `SigGenSetup` object from `MJDSigGen.jl`.

# The same function as in the SSD example can be used to obtain fieldgen simulated capacitance matrix

LegendGeSim.capacitance_matrix(sim_fieldgen)

# ### Alternative

# Alternative input using an already existing `impurity_profile` and `offset_in_mm` as input fields

simulation_settings_siggen1 = Dict(
    "method" => "fieldgen",
    "cached_name" => "", # a non-empty string will cache the simulation result
    "impurity_profile" => "cache/V99000.dat",
    "offset_in_mm" => 0
);

#

sim_fieldgen = LegendGeSim.simulate_fields(detector_metadata_filename, environment_settings, simulation_settings_siggen1; overwrite=true)

# ## Compare `SSD` to `Fieldgen`

# We can construct wlectric and weighting potential based on Fieldgen results in the same format as SSD

e_pot = SSD.ElectricPotential(sim_fieldgen);
w_pot = SSD.WeightingPotential(sim_fieldgen);

# ToDo: same units for SSD and fieldgen
plot(
    plot(sim_ssd.electric_potential, full_det = true, title = "Electric potential via SSD"),
    plot(e_pot, full_det = true, title = "Electric potential via Fieldgen"),
    size = (1000, 500), layout = (1, 2)
)

#jl savefig("SSD_vs_fieldgen_Epot.pdf") # hide
#md savefig("SSD_vs_fieldgen_Epot.svg"); nothing # hide
#md # ![SSD_vs_fieldgen_Epot](SSD_vs_fieldgen_Epot.svg)

plot(
    plot(sim_ssd.weighting_potentials[1], full_det = true, title = "Weighting potential via SSD"),
    plot(w_pot, full_det = true, title = "Weighting potential via fieldgen"),
    size = (1000, 500), layout = (1, 2)
)

#jl savefig("SSD_vs_fieldgen_Wpot.pdf") # hide
#md savefig("SSD_vs_fieldgen_Wpot.svg"); nothing # hide
#md # ![SSD_vs_fieldgen_Wpot](SSD_vs_fieldgen_Wpot.svg)
