"""
    simulate_detector(sim_config_filename)

AbstractString -> <detector simulation>

Simulate detector according to given simulation configuration.
What type <detector simulation> is depends on the given method of simulation.
"""
function simulate_detector(sim_config::PropDict)
    det_meta = propdict(sim_config.detector_metadata)
    env = Environment(sim_config)
    ps_simulator = PSSimulator(sim_config)

    simulate_detector(det_meta, env, sim_config.pss.cached_name, ps_simulator)
end


function simulate_detector(sim_config_filename::AbstractString)
    simulate_detector(propdict(sim_config_filename))
end


function simulate_detector(det_metadata::AbstractString, sim_config_filename::AbstractString)
    # merge detector info with the rest 
    sim_config = load_config(det_metadata, sim_config_filename)
    simulate_detector(sim_config)
end


"""
    simulate_detector(det_meta, env, config_name, ::SSDSimulator)

PropDict, Environment, AbstractString -> SSD.Simulation 

Look up cached SSD simulation .h5f file, and if does not exist,
    simulate the detector according to geometry as given in LEGEND metadata <det_meta>
    and environment settings given in <env> using SolidStateDetectors.
The simulation will be cached based on <config_name>.

"""
function simulate_detector(det_meta::PropDict, env::Environment, cached_name::AbstractString, ssd_sim::SSDSimulator)
    det_h5 = joinpath("cache", cached_name*"_ssd.h5f")
    if isfile(det_h5)
        @info("Reading SSD simulation from cached file $det_h5")
        simulation = SolidStateDetectors.ssd_read(det_h5, Simulation)
    else
        @info("Simulating $cached_name with SSD from scratch for given settings")
        ssd_conf = ssd_config(det_meta, env, ssd_sim)
        simulation = simulate_ssd(ssd_conf)

        if !ispath(dirname(det_h5)) mkpath(dirname(det_h5)) end
        SolidStateDetectors.ssd_write(det_h5, simulation)
        @info("-> Saved cached simulation to $det_h5")
    end

    simulation
end    


"""
    simulate_detector(det_meta, env, config_name, ps_simulator)
    
PropDict, Environment, AbstractString, Siggen -> SigGenSetup       

Construct a fieldgen/siggen configuration file
    according to geometry as given in LEGEND metadata <det_meta>
    and environment settings given in <env>.
Look up fieldgen generated electric field and weighting potential files
    based on given <config_name>, and if not found, call fieldgen.
"""
function simulate_detector(det_meta::PropDict, env::Environment, cached_name::AbstractString, ps_simulator::SiggenSimulator)
    @info "Constructing Fieldgen/Siggen configuration file"
    # returns the name of the resulting siggen config file
    # and the name of (already or to be) generated weighting potential file
    siggen_config_name, fieldgen_wp_name = siggen_config(det_meta, env, ps_simulator, cached_name)
    fieldgen_wp_name = joinpath("cache", fieldgen_wp_name)

    # if the WP file with such a name exists...
    if isfile(fieldgen_wp_name)
        #...do nothing, siggen will later read the files based on the generated conifg
        @info "Reading cached fields from $fieldgen_wp_name"
    else
        #...otherwise call fieldgen
        @info "_|~|_|~|_|~|_ Fieldgen simulation"
        fieldgen(siggen_config_name)
        @info "_|~|_|~|_|~|_ Fieldgen simulation complete"
    end

    # a SigGenSetup object
    SigGenSetup(siggen_config_name)            
end

# ------------------------------------------------------------------- SSD

"""
    simulate_ssd(ssd_config)

Dict -> SSD.Simulation 

Simulate detector based on given SSD config
"""
function simulate_ssd(ssd_config::Dict)

    @info "_|~|_|~|_|~|_ SSD simulation"

    simulation = Simulation(SolidStateDetector{Float32}(ssd_config))

    println("...electric potential")
    calculate_electric_potential!( simulation,
                               max_refinements = 4)

    println("...electric field")
    calculate_electric_field!(simulation, n_points_in_φ = 72)

    println("...drift field")
    calculate_drift_fields!(simulation)

    println("...weighting potential")
    for contact in simulation.detector.contacts
        calculate_weighting_potential!(simulation, contact.id, max_refinements = 4, n_points_in_φ = 2, verbose = false)
    end

    @info "_|~|_|~|_|~|_ SSD detector simulation complete"

    simulation
end


"""
    ssd_config(meta, env)

PropDict, Environment -> Dict 

Construct SSD type config based on
    geometry as given in LEGEND metadata <meta>
    and environment settings given in <env>.
""" 
function ssd_config(meta::PropDict, env::Environment, ssd_sim::SSDSimulator)
    ## some parameters that will be used multiple times later
    # or need a calculation
    cylinder_height = meta.geometry.height_in_mm
    cylinder_radius = meta.geometry.radius_in_mm
    borehole_height = cylinder_height - meta.geometry.well.gap_in_mm
    borehole_r_bottom = meta.geometry.well.radius_in_mm
    borehole_r_top = borehole_r_bottom + borehole_height * tan(meta.geometry.taper.top.inner.angle_in_deg / 180. * π)
    cone_r_top = cylinder_radius - meta.geometry.taper.top.outer.height_in_mm * tan(meta.geometry.taper.top.outer.angle_in_deg / 180. * π)
    cone_height = meta.geometry.taper.top.outer.height_in_mm

    ## grid 
    # kind of initialize so that it exists out of the loop?    
    grid = Dict()

    if ssd_sim.coord == "cylindrical"
        grid = Dict(
            "axes" => Dict(
                "r" => Dict("to" => cylinder_radius*1.2, "boundaries" => "inf"), # leave some margin
                "phi" => Dict("from" => 0, "to" => ssd_sim.comp == "3D" ? 360 : 0, "boundaries" => "periodic"),
                "z" => Dict(
                    "from" => -0.2*cylinder_height,
                    "to" => cylinder_height*1.2, # leave some margin
                    "boundaries" => Dict("left" => "inf", "right" => "inf")
                )
            )
        )
    else # cartesian
        grid = Dict(
            "axes" => Dict(
                "x" => Dict("from" => -cylinder_radius*1.2, "to" => cylinder_radius*1.2, "boundaries" => "inf"),
                "y" => Dict("from" => -cylinder_radius*1.2, "to" => cylinder_radius*1.2, "boundaries" => "inf"),
                "z" => Dict("from" => -0.2*cylinder_height, "to" => cylinder_height*1.2, "boundaries" => "inf")
            )
        )
    end

    grid = merge(grid, Dict("coordinates" => ssd_sim.coord))

    ## charge density model
    # kind of initialize so that it exists out of the loop?    
    charge_density_model = Dict()

    if ssd_sim.coord == "cylindrical"
        charge_density_model = Dict(
            "r" => Dict("init" => 0, "gradient" => 0),
            "phi" => Dict("init" => 0.0, "gradient" => 0.0),
            "z" => Dict("init" => -1e7, "gradient" => 1e5)
        )
    else # cartesian
        charge_density_model = Dict(
            "x" => Dict("init" => 0.0, "gradient" => 0.000),
            "y" => Dict("init" => 0.0, "gradient" => 0.0),
            "z" => Dict("init" => 1e7, "gradient" => 5e4)
        )
    end

    # for now we only use linear charge density model for either cartesian or cylindrical coordinates
    # charge_density_model["name"] = "linear"
    charge_density_model = merge(charge_density_model, Dict("name" => "linear"))

    ## geometry
    dct = Dict(
        "name" => meta.det_name,
        "medium" => env.medium,
        "units" => Dict(
            "length" => "mm",
            "angle" => "deg",
            "potential" => "V",
            "temperature" => "K"
        ),
        "grid" => grid,
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
                        # crystal
                        Dict(
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
                                    "h" => cylinder_height,
                                    "translate" => Dict("z" => 0.0)
                                )#, # cylinder
                                # upper cone
                                Dict(
                                    "name" => "Upper cone",
                                    "type" => "cone",
                                    "r" => Dict(
                                        "bottom" => Dict("from" => cylinder_radius, "to" => cylinder_radius+1),
                                        "top" => Dict("from" => cone_r_top, "to" => cylinder_radius+1)
                                    ),
                                    "phi" => Dict("from" => 0.0, "to" => 360.0),
                                    "h" => meta.geometry.taper.top.outer.height_in_mm,
                                    "translate" => Dict("z" => cylinder_height - cone_height)
                                ) # upper cone
                            ] # crystal parts
                        ), # crystal
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
                    "parts" => [Dict() for i in 1:6] # to be filled later
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
        "h" => cylinder_height - cone_height
    )

    # upper cone 
    dct["objects"][3]["geometry"]["parts"][3] = Dict(
        "name" => "Upper cone",
        "type"=> "cone",
        "r" => Dict(
            "bottom" => Dict("from" => cylinder_radius, "to" => cylinder_radius),
            "top" => Dict("from" => cone_r_top, "to" => cone_r_top)
        ),
        "phi" => Dict("from" => 0.0, "to" => 360.0),
        "h" => cone_height,
        "translate" => Dict("z" => cylinder_height - cone_height)
    )

    # top lid
    dct["objects"][3]["geometry"]["parts"][4] = Dict(
        "name" => "top lid",
        "type" => "tube",
        "r" => Dict(
            "from" => borehole_r_top, 
            "to" => cone_r_top
        ),
        "phi" => Dict("from" => 0.0, "to" => 360.0),
        "h" => 0.0,
        "translate" => Dict("z" => cylinder_height)
    )

    # borehole lining
    dct["objects"][3]["geometry"]["parts"][5] = Dict(
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
    dct["objects"][3]["geometry"]["parts"][6] = Dict(
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


"""
    ssd_config(det_meta)

AbstractString -> Dict 

Construct SSD type config based on geometry as given in LEGEND metadata <det_meta>.
Set arbitrary envorinment variables and simulation settings (used only for geometry visualization)
"""
function ssd_config(det_meta::AbstractString)
    ssd_config(propdict(det_meta), Environment(), SSDSimulator())
end


"""
    siggen_config(meta, env, siggen_sim, config_name)

PropDict, Environment, SiggenSimulator, String -> String, String    

Construct Siggen/Fieldgen config based on
    geometry as given in LEGEND metadata <meta>,
    environment settings given in <env>,
    and fieldgen settings contained in <siggen_sim>.

The resulting siggen config will be cached based on given <config_name>
The function returns the path to the cached siggen config file, and
    the relative path to the fieldgen WP file (to check if already exists later)
"""
function siggen_config(meta::PropDict, env::Environment, siggen_sim::SiggenSimulator, cached_name::AbstractString)
    # function pss_config(meta::PropDict, env::Environment, siggen_sim::SiggenSimulatorss, config_name::AbstractString)
    # construct siggen geometry based on LEGEND metadata JSON + extra environment settings
    println("...detector geometry")
    siggen_geom = meta2siggen(meta, env)

    # construct lines for future siggen config file 
    detector_lines = Vector{String}()
    push!(detector_lines, "# detector related parameters")
    for (param, value) in siggen_geom
        push!(detector_lines, "$param   $value")
    end

    println("...fieldgen configuration")
    # create filenames for future/existing fieldgen output 
    fieldgen_wp_name = joinpath("fieldgen", cached_name*"_fieldgen_WP.dat")
    fieldgen_names = Dict(
        "drift_name" => joinpath("..", siggen_sim.drift_vel),
        # "drift_name" => joinpath("..", "data", "drift_vel_tcorr.tab"),
        "field_name" => joinpath("fieldgen", cached_name*"_fieldgen_EV.dat"),
        "wp_name" => fieldgen_wp_name
    )

    # add corresponding lines to existing detector geometry lines
    push!(detector_lines, "")
    push!(detector_lines, "# fieldgen filenames")
    for (param, value) in fieldgen_names
        push!(detector_lines, "$param   $value")
    end
    
    # read lines from user input for fieldgen
    fieldgen_lines = readlines(open(siggen_sim.fieldgen_config))
    fieldgen_lines = replace.(fieldgen_lines, "\t" => "   ")

    # unite constructed detector geometry and fieldgen settings
    total_lines = vcat(detector_lines, [""], fieldgen_lines)

    # write a siggen/fieldgen config file 
    sig_conf_name = joinpath("cache", cached_name*"_siggen_config.txt")
    writedlm(sig_conf_name, total_lines)
    println("...fieldgen/siggen config written to $sig_conf_name")

    # return resulting fieldgen/siggen config name 
    # (currently kinda redundant but that's because we have no idea what's gonna happen later)
    sig_conf_name, fieldgen_wp_name
end


"""
    meta2siggen(meta, env)

PropDict, Environment -> Dict    

Construct siggen type config in a form of Dict based on
    geometry as given in LEGEND metadata <meta>
    and environment settings given in <env>.
"""
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

    # comment from David Radford:

    # The bottom taper is always at 45 degrees, but the top outer taper is not.
    # It's width/angle is defined using either the outer_taper_width parameter or the taper_angle parameter.
    # Likewise, the inner taper width/angle is defined using either the inner_taper_width parameter or the taper_angle parameter.
    # The inner_taper_length can be anything from zero to the hole_length
    # TODO: check if bottom taper angle is 45 deg, if not issue a warning for siggen case

        # size of 45-degree taper at bottom of ORTEC-type crystal (for r=z)
        "bottom_taper_length" => meta.geometry.taper.bottom.outer.height_in_mm,
        # z-length of outside taper for inverted-coax style
        "outer_taper_length" => meta.geometry.taper.top.outer.height_in_mm,
        # z-length of inside (hole) taper for inverted-coax style
        "inner_taper_length" => meta.geometry.taper.top.inner.height_in_mm,
        # taper angle in degrees, for inner or outer taper
        "taper_angle" => meta.geometry.taper.top.outer.angle_in_deg,
        # "taper_angle" => meta.geometry.taper.top.inner.angle_in_deg,
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




