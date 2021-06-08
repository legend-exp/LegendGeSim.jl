function detector_config(det_metadata::AbstractString, env::Environment, ps_simulator::PSSimulator)
    meta = PropDicts.read(PropDict, det_metadata)
    detector_config(meta, env, ps_simulator)
end    


"""
    SSDconfig(det_metadata)

Read LEGEND metadata json file and construct a dictionary
with SSD detector configuration format
"""
function detector_config(meta::PropDict, env::Environment, ::SSDSimulator)
    cylinder_height = meta.geometry.height_in_mm
    cylinder_radius = meta.geometry.radius_in_mm
    borehole_height = cylinder_height - meta.geometry.well.gap_in_mm
    borehole_r_bottom = meta.geometry.well.radius_in_mm
    borehole_r_top = borehole_r_bottom + borehole_height * tan(meta.geometry.taper.top.inner.angle_in_deg / 180. * π)

    # I'm not sure if this setting is common to all detectors
    # Currently working with V02160A
    charge_density_model = Dict(
        "name" => "linear",
        "r" => Dict(
            "init" => 0,
            "gradient" => 0
        ),
        "phi" => Dict(
            "init" => 0.0,
            "gradient" => 0.0
        ),
        "z" => Dict(
            "init" => -1e7,
            "gradient" => 1e5
        )
    )

    dct = Dict(
        "name" => meta.det_name,
        "units" => Dict(
            "length" => "mm",
            "angle" => "deg",
            "potential" => "V",
            "temperature" => "K"
        ),
        "objects" => [
            # crystal
            Dict(
                "type" => "semiconductor",
                "material" => "HPGe",
                "bulk_type" => "p",
                "temperature" => env.crystal_t,
                "charge_density_model" => charge_density_model,
                "geometry" => Dict(
                    "type" => "difference",
                    "parts" => [
                        # cylinder
                        Dict(
                            "name" => "Initial Cylinder",
                            "type" => "tube",
                            "r" => Dict(
                                "from" => 0.0,
                                "to" => cylinder_radius,
                            ),
                            "phi" => Dict("from" => 0.0, "to" => 360.0),
                            "h" => cylinder_height
                        ), # cylinder
                        # borehole
                        Dict(
                            "name" => "borehole",
                            "type" => "cone",
                            "r" => Dict(
                                "bottom" => Dict(
                                    "from" => 0.0,
                                    "to" => borehole_r_bottom
                                ),
                                "top" => Dict(
                                    "from" => 0.0,
                                    "to" => borehole_r_top
                                )
                            ),
                            "phi" => Dict("from" => 0.0, "to" => 360.0),
                            "h" => borehole_height,
                            "translate" => Dict("z" => meta.geometry.well.gap_in_mm)
                        ), # borehole
                        # groove
                        Dict(
                            "name" => "ditch",
                            "type" => "tube",
                            "r" => Dict(
                                "from" => meta.geometry.groove.outer_radius_in_mm - meta.geometry.groove.width_in_mm,
                                "to" => meta.geometry.groove.outer_radius_in_mm
                            ),
                            "phi" => Dict("from" => 0.0, "to" => 360.0),
                            "h" => meta.geometry.groove.depth_in_mm
                        ) # groove
                    ] # parts
                ) # geometry
            ), # crystal
            # p+ contact 
            Dict(
                "type" => "contact",
                "material" => "HPGe",
                "channel" => 1,
                "potential" => 0.0,
                "geometry" => Dict(
                    "type" => "tube",
                    "r" => Dict(
                        "from" => 0.0,
                        "to" => meta.geometry.contact.radius_in_mm
                    ),
                    "phi" => Dict("from" => 0.0, "to" => 360.0),
                    "h" => meta.geometry.contact.depth_in_mm
                )
            ), # p+ contact
            # n+ contact
            Dict(
                "type" => "contact",
                "material" => "HPGe",
                "channel" => 2,
                "potential" => env.op_voltage == 0 ? meta.characterization.manufacturer.op_voltage_in_V : env.op_voltage,
                "geometry" => Dict(
                    "type" => "union",
                    "parts" => [Dict() for i in 1:5] # to be filled later
                )
            ) # n+ contact
        ] # objects
    )

    ## Construct n+ surface

    # bottom lid
    dct["objects"][3]["geometry"]["parts"][1] = Dict(
        "name" => "bottom lid",
        "type" => "tube",
        "r" => Dict(
            "from" => meta.geometry.groove.outer_radius_in_mm,
            "to" => cylinder_radius # cylinder radius
        ),
        "phi" => Dict("from" => 0.0, "to" => 360.0),
        "h" => 0
    )

    # side lining
    dct["objects"][3]["geometry"]["parts"][2] = Dict(
        "name" => "side lining",
        "type" => "tube",
        "r" => Dict(
            "from" => cylinder_radius, 
            "to" => cylinder_radius
        ),
        "phi" => Dict("from" => 0.0, "to" => 360.0),
        "h" => cylinder_height
    )

    # top lid
    dct["objects"][3]["geometry"]["parts"][3] = Dict(
        "name" => "top lid",
        "type" => "tube",
        "r" => Dict(
            "from" => borehole_r_top, 
            "to" => cylinder_radius 
        ),
        "phi" => Dict("from" => 0.0, "to" => 360.0),
        "h" => 0.0,
        "translate" => Dict("z" => cylinder_height)
    )

    # borehole lining
    dct["objects"][3]["geometry"]["parts"][4] = Dict(
        "name" => "borehole linning",
        "type" => "cone",
        "r" => Dict(
            "bottom" => Dict(
                "from" => borehole_r_bottom, 
                "to" => borehole_r_bottom
            ),
            "top" => Dict(
                "from" => borehole_r_top,
                "to" => borehole_r_top
            )
        ),
        "phi" => Dict("from" => 0.0, "to" => 360.0),
        "h" => borehole_height,
        "translate" => Dict("z" => meta.geometry.well.gap_in_mm)
    )

    # borehole bottom lid
    dct["objects"][3]["geometry"]["parts"][5] = Dict(
        "name" => "borehole bottom lid",
        "type" => "tube",
        "r" => Dict(
            "from" => 0.0,
            "to" => borehole_r_bottom
        ),
        "phi" => Dict("from" => 0.0, "to" => 360.0),
        "h" => 0.0,
        "translate" => Dict("z" => meta.geometry.well.gap_in_mm)
    )

    dct
end

# """
# Create SiggenSetup object based on given full siggen config
# (geometry + fieldgen)
# """
# function detector_config(::PropDict, ::Environment, siggen_sim::SiggenFull)
#     SigGenSetup(siggen_sim.siggen_config)    
# end


"""
Construct siggen config file based on detector geometry given via LEGEND metadata JSON
and given fieldgen settings (contained within SiggenFieldgen object)
"""
function detector_config(meta::PropDict, env::Environment, siggen_sim::SiggenFieldgen)
    # construct siggen geometry based on LEGEND metadata JSON + extra environment settings
    siggen_geom = meta2siggen(meta, env)

    # construct lines for siggen config file 
    detector_lines = Vector{String}()
    for (param, value) in siggen_geom
        push!(detector_lines, "$param   $value")
    end


    # read lines from user input for fieldgen 
    fieldgen_lines = readlines(open(siggen_sim.fieldgen_config))
    fieldgen_lines = replace.(fieldgen_lines, "\t" => "  ")

    total_lines = vcat(["# detector related parameters"], detector_lines, fieldgen_lines)

    # write a config file 
    # not possible because siggen is looking in the same path as the config, which becomes /tmp instead of data/
    # sig_conf_name = mktemp()[1]
    sig_conf_name = siggen_sim.fieldgen_config * ".temp"
    writedlm(sig_conf_name, total_lines)

    # construct setup
    SigGenSetup(sig_conf_name)    
end


function meta2siggen(meta::PropDict, env::Environment)
    Dict(
    # crystal
        # z length: crystal height
        "xtal_length" => meta.geometry.height_in_mm,
        # radius: crystal radius
        "xtal_radius" => meta.geometry.radius_in_mm,
    # contact
        # point contact length
        "pc_length" => meta.geometry.contact.depth_in_mm,
        # point contact radius
        "pc_radius" => meta.geometry.contact.radius_in_mm,
    # groove
        # wrap-around radius. Set to zero for ORTEC: groove outer radius
        "wrap_around_radius" => meta.geometry.groove.outer_radius_in_mm,
        # depth of ditch next to wrap-around for BEGes. Set to zero for ORTEC: groove height
        "ditch_depth" => meta.geometry.groove.depth_in_mm,
        # width of ditch next to wrap-around for BEGes. Set to zero for ORTEC: groove inner radius
        "ditch_thickness" => meta.geometry.groove.width_in_mm,
    # borehole 
        # length of hole, for inverted-coax style
        "hole_length" => meta.geometry.height_in_mm - meta.geometry.well.gap_in_mm,
        # radius of hole, for inverted-coax style
        "hole_radius" => meta.geometry.well.radius_in_mm,
    # tapering
        # size of 45-degree taper at bottom of ORTEC-type crystal (for r=z)
        "bottom_taper_length" => meta.geometry.taper.bottom.outer.height_in_mm,
        # z-length of outside taper for inverted-coax style
        "outer_taper_length" => meta.geometry.taper.top.outer.height_in_mm,
        # z-length of inside (hole) taper for inverted-coax style
        "inner_taper_length" => meta.geometry.taper.top.inner.height_in_mm,
        # taper angle in degrees, for inner or outer taper
        "taper_angle" => meta.geometry.taper.top.inner.angle_in_deg, # still confused whether this angle means outer or inner taper,
    # depth of full-charge-collection boundary for Li contact (not currently used)
        "Li_thickness" => meta.geometry.dl_thickness_in_mm,

    # configuration for mjd_fieldgen (calculates electric fields & weighing potentials)
        # detector bias for fieldgen, in Volts
        "xtal_HV" => env.op_voltage == 0 ? meta.characterization.manufacturer.op_voltage_in_V : env.op_voltage,

    # configuration for signal calculation 
        # crystal temperature in Kelvin
        "xtal_temp" => env.crystal_t
    )
end


# function siggen2meta(siggen_config::AbstractString)

#     # read lines from file
#     siggen_lines = readlines(open(siggen_config))

#     # construct dict
#     siggen_dct = Dict()
#     for line in siggen_lines
#         if not "=" in line continue end
#         param, value = split(line, "=")
#         siggen_dct[strip(param)] = parse(Float32, strip(value))
#     end

#     siggen2meta(siggen_dct)
# end


# function siggen2meta(siggen_dct::Dict)

#     ## Construct LEGEND metadata Dict

#     dct_meta["geometry"] = Dict(
#         "contact" => Dict(),
#         "groove" => Dict(),
#         "well" => Dict(),
#         "taper" => Dict(
#             "bottom" => Dict("outer" => Dict()),
#             "top" => Dict(
#                 "outer" => Dict(),
#                 "inner" => Dict()
#             )
#         )
#     )

#     # z length: crystal height
#     dct_meta["geometry"]["height_in_mm"] = siggen_dct["xtal_length"]
#     # radius: crystal radius
#     dct_meta["geometry"]["radius_in_mm"] = siggen_dct["xtal_radius"]

#     # contact
#     # point contact length
#     dct_meta["geometry"]["contact"]["depth_in_mm"] = siggen_dct["pc_length"]
#     # point contact radius
#     dct_meta["geometry"]["contact"]["radius_in_mm"] = siggen_dct["pc_radius"]

#     # groove
#     # wrap-around radius. Set to zero for ORTEC: groove outer radius
#     dct_meta["geometry"]["groove"]["outer_radius_in_mm"] = siggen_dct["wrap_around_radius"]
#     # depth of ditch next to wrap-around for BEGes. Set to zero for ORTEC: groove height
#     dct_meta["geometry"]["groove"]["depth_in_mm"] = siggen_dct["ditch_depth"]
#     # width of ditch next to wrap-around for BEGes. Set to zero for ORTEC: groove inner radius
#     dct_meta["geometry"]["groove"]["width_in_mm"] = siggen_dct["ditch_thickness"]
    
#     # borehole 
#     # length of hole, for inverted-coax style
#     dct_meta["geometry"]["well"]["gap_in_mm"] = siggen_dct["hole_length"]
#     # radius of hole, for inverted-coax style
#     dct_meta["geometry"]["well"]["radius_in_mm"] = siggen_dct["hole_length"]
    
#     # tapering
#     # size of 45-degree taper at bottom of ORTEC-type crystal (for r=z)
#     dct_meta["geometry"]["taper"]["bottom"]["outer"]["height_in_mm"] = siggen_dct["bottom_taper_length"]
#     # z-length of outside taper for inverted-coax style
#     dct_meta["geometry"]["taper"]["top"]["outer"]["height_in_mm"] = siggen_dct["outer_taper_length"]
#     # z-length of inside (hole) taper for inverted-coax style
#     dct_meta["geometry"]["taper"]["top"]["inner"]["height_in_mm"] = siggen_dct["inner_taper_length"]
#     # taper angle in degrees, for inner or outer taper
#     dct_meta["geometry"]["taper"]["top"]["inner"]["angle_in_mm"] = siggen_dct["taper_angle"]

#     # depth of full-charge-collection boundary for Li contact (not currently used)
#     dct_meta["geometry"]["dl_thickness_in_mm"] = siggen_dct["Li_thickness"]

#     ## Construct Environment struct

#     env = Environment(siggen_dct["xtal_temp"], siggen_dct["xtal_HV"])

#     PropDict(dct_meta), env

# end



"""
    simulate_detector(det_path, det_name)

Simulate detector based on json geometry.

det_path: path to detector json file
det_name: detector name (e.g. "V05266A").
The code will look for a .json file "det_path/det_name.json"

Output: SSD detector simulation object
"""
function simulate_detector(detector_config::Dict)

    simulation = Simulation(SolidStateDetector{Float32}(detector_config))

    @info("-> Electric potential...")
    calculate_electric_potential!( simulation,
                               max_refinements = 4)

    @info("-> Electric field...")
    calculate_electric_field!(simulation, n_points_in_φ = 72)

    @info("-> Capacitance...")
    calculate_capacitance(simulation)

    @info("-> Drift field...")
    calculate_drift_fields!(simulation)

    @info("-> Weighting potential...")
    for contact in simulation.detector.contacts
        calculate_weighting_potential!(simulation, contact.id, max_refinements = 4, n_points_in_φ = 2, verbose = false)
    end

    @info "Detector simulation complete"

    det_h5 = joinpath("cache", basename(detector_config["name"]*".h5f"))

    if !ispath(dirname(det_h5)) mkpath(dirname(det_h5)) end
    SolidStateDetectors.ssd_write(det_h5, simulation)

    @info("-> Saved cached simulation to $det_h5")

    simulation
end