# Field Simulation

As inputs to field simulation you should provide

1. LEGEND detector metadata

2. Environment settings

3. Simulation settings

## Example

```julia
detector_metadata_filename = "path/to/legend-detectors/germanium/diodes/V04545A.json"

environment_settings = Dict(
    "crystal_temperature_in_K" => 77, 
    "medium" => "vacuum", 
    "dl" => "vendor" # optional, default 0
    "operating_voltage_in_V" => 5000, # optional, default recV from metadata
)

simulation_settings = Dict(
    "method" => "SSD",
    "cached_name" => "test",
    "crystal_metadata_path" => "path/to/legend-detectors/germanium/crystals"
)

sim = LegendGeSim.simulate_fields(detector_metadata_filename, environment_settings, simulation_settings; overwrite=true)
```

Provide an optional argument `overwrite` to overwrite an existing file with the same `cached_name`. Default: `false`

In case of SSD as method, the output of `simulate_fields()` is a `Simulation` object from `SolidSateDetectors.jl`. Among other things, it contains a `SolidStateDetector` object under `sim.detector`. See the short tutorial on how to [Visualize Detector Geometry](@ref) using this object.

In case of siggen as method, the output of `simulate_fields()` is a `SigGenSetup` object from `MJDSigGen.jl` containing information needed for the siggen simulation. It is possible to extract the electric and weighting potential from this object and plot it (see tutotorial [Simulate Fields](@ref)).

## 1. LEGEND detector metadata

To obtain JSONs for all LEGEND detectors you can clone the repository `legend-detectors`. The JSONs can be found under `germanium/diodes`. You may also create your custom detector JSON file as long as it follows the LEGEND detector metadata format (see README in `diodes/`).

## 2. Environment settings

Environment settings contain:
- `"crystal_temperature_in_K"`: crystal temperature in units of Kelvin e.g. 77 for liquid nitrogen, 90 for LAr
- `"medium"`: `"vacuum"` or `"LAr"` (note: the LAr option has actually never been tested yet)
- `"operating_voltage_in_V"`: bias voltage at the n+ contact in units of Volts

    If you do not provide this field in the environment settings, the value from `characterization.l200_site.recommended_voltage_in_V` in detector metadata JSON will be taken.

- `"dl"`: dead layer (DL) thickness in mm.

    This is not really part of "environment" but is added as a temporary quickfix for now to be able to study DL effects. In the future might be moved to simulation settings (see below).

    Provide a value e.g. `"dl" => 1.2`. If you do not provide this setting at all, default `0` will be taken
    
    Providing `"dl" => "vendor"` means the value from `characterization.manufacturer.dl_thickness_in_mm` in detector metadata JSON file will be taken.

## 3. Simulation settings

Simulation settings contain:

- `"method"`: `"SSD"` (or `"ssd"`) or `"siggen"` (or `"fieldgen"`)

    Based on the given method, `LegendGeSim` will automatically construct configurations for either method based on the detector metadata, environment, and simulation settings

- `"cached_name"`: if provided, your simulation will be cached in a file for future uses in a directory `cache/` (will be created if not present)

    Format:

    - `cache/<detector name>_<cached_name>_SSD.h5f` file for SSD

    - `cache/<detector_name>_<cached_name>_fieldgen_WP.dat` and `*_EV.dat` files for siggen.
        
        In case of siggen, a siggen config named `cache/<detector_name>_<cached_name>_siggen_config.txt` will also be saved (generated based on given settings and geometry).
        
    The detector name is taken from the detector metadata JSON `"name"` field.
        
    If cached name is not provided (or `""` is provided), the simulation will not be cached (for siggen will create a temporary config file).
        
    It may be convenient to compose the `cached_name` in the simulation settings based on environment settings, varying detector geometry etc.

- `"crystal_metadata_path"`: path to folder with crystal metadata JSONs.

    To obtain JSONs for LEGEND crystals you can clone the repository `legend-detectors`. The JSONs can be found under `germanium/crystals`. You may also create your custom crystal JSON file as long as it follows the LEGEND crystal metadata format (see README in `crystals/`).

    The crystal corresponding to your detector will be chosen from the given path based on the `"name"` field in the detector JSON. E.g. if your detector `"name"` field states `"V04545A"`, `LegendGeSim` will be looking for a file `"V04545.json"`.
    
    The impurity profile will be constructed by fitting the impurity measurement points in the JSON using David Radford's function `a + b*z + c*exp((L - z)/tau)` (in crystal axis coordinates, later converted to detector axis). In case of SSD, the profile will be provided to the simulation internally, while for siggen, an input file will be created in siggen format (e.g. `"cache/V04545.dat"`). In either case, the position of the detector in the crystal is determined by the offset field contained in the crystal metadata.

    For SSD, if no crystal metadata path is provided (or `""` is provided), a dummy constant impurity of 10⁹ e/cm³ will be used.

Alternative options for `siggen`: rather than providing a path to the crystal metadata folder and implementing impurity profile as described above, you may instead provide two inputs as is traditionally done when running siggen:

- `"impurity_profile"`: path to a `.dat` file with the impurity profile (unformatted stream of Float32 of impurity values in units of 10¹⁰ e/cm³ with a step of 1 mm from tail end)

- `"offset_in_mm"`: offset of the detector from tail end in mm

If neither `"crystal_metadata_path"` not these settings are provided for siggen, a dummy constant impurity will be used as described above.

Optional settings for `siggen`:

- `"fieldgen_config"`: path to a siggen format config with "extra" settings e.g. crystal grid. Default `examples/configs/fieldgen_settings.txt` (in `LegendGeSim`). Note: in case of no crystal impurity profile input, a dummy constant impurity that is used is written in this file. In the future, all other optional fieldgen settings will be incorporated into the simulation settings input via `LegendGeSim`.

- `"drift_vel"`: path to siggen format drift velocity file. Default `examples/configs/drift_vel_tcorr.tab` (in `LegendGeSim`)

**WARNING**

The implementation here is kind of janky, and if you want to read from a cached file, in theory you should only need to provide the cached name. However, for now it will also ask you to provide the other non-optional fields as well, and even worse - it will print out those settings first, and only then figure out that there's a cached file.

This will of course be fixed in the future to not require those inputs, and possibly issue a warning if a cache file exists but settings are provided etc.

However, even after that, if you provide certain inputs like temperature and crystal metadata path and a cached name that already exists, that cached file will be read even if the simulation there does no correspond to your other inputs, unless you provide `overwrite=true` to overwrite that file. For now at least, it is the user's job to remember which settings were used for the cached simulation.

**User suggestions are very welcome on how to manage this better**    