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


function detector_config(det_metadata::AbstractString, ps_simulator::PSSimulator)
    meta = PropDicts.read(PropDict, det_metadata)
    detector_config(meta, ps_simulator)
end    


"""
    SSDconfig(det_metadata)

Read LEGEND metadata json file and construct a dictionary
with SSD detector configuration format
"""
function detector_config(meta::PropDict, ssd_sim::SSDSimulator)
    cylinder_height = meta.geometry.height_in_mm
    cylinder_radius = meta.geometry.radius_in_mm
    borehole_height = cylinder_height - meta.geometry.well.gap_in_mm
    borehole_r_bottom = meta.geometry.well.radius_in_mm
    borehole_r_top = borehole_r_bottom + borehole_height * tan(meta.geometry.taper.top.inner.angle_in_deg / 180. * π)

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
                "temperature" => ssd_sim.crystal_t,
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
                "potential" => meta.characterization.manufacturer.op_voltage_in_V,
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


# function detector_config(meta::PropDict, siggen_sim::SiggenSimulator)
    # # setup object with settings that are not in metadata
    # setup = SigGenSetup()
    # # setup = SigGenSetup(siggen_sim.fieldgen_input)

    # setup.verbosity = 0   #  0 = terse, 1 = normal, 2 = chatty/verbose

    # # configuration for mjd_fieldgen (calculates electric fields & weighing potentials)
    # setup.xtal_grid = 0.1		# grid size in mm for field files (usually 0.5 or 0.1 mm)
    # setup.impurity_z0 = -0.92	# net impurity concentration at Z=0, in 1e10 e/cm3
    # setup.impurity_gradient = 0.031	# net impurity gradient, in 1e10 e/cm4
    # setup.impurity_quadratic = 0.01	# net impurity difference from linear, at z=L/2, in 1e10 e/cm3

    # # options for mjd_fieldgen:
    # setup.max_iterations = 30000  # maximum number of iterations to use in mjd_fieldgen
    # setup.write_field = 1      # 0/1: do_not/do write the standard field output file
    # setup.write_WP = 1      # 0/1: do_not/do calculate the weighting potential and write it to the file

    # # file names
    # setup.drift_name = "drift_vel_tcorr.tab"    # drift velocity lookup table
    # setup.field_name = "fields/ev_GWD6022.dat"    # potential/efield file name; no included spaces allowed
    # setup.wp_name = "fields/wp_GWD6022.dat"    # weighting potential file name; no included spaces allowed

    # # configuration for signal calculation
    # setup.preamp_tau = 0		# integration time constant for preamplifier, in ns
    # setup.time_steps_calc = 40000	# max. number of time steps used in calculations
    # setup.step_time_calc = 0.1	# length of time step used for calculation, in ns
    # setup.step_time_out = 1.0	# length of time step for output signal, in ns

    # # nonzero values in the next few lines significantly slows down the code
    # setup.energy = 0
    # setup.charge_cloud_size = 0      # initial FWHM of charge cloud, in mm
    # setup.use_diffusion = 0      # set to 0/1 for ignore/add diffusion as the charges drift
    # setup.use_acceleration = 0
    # setup.use_repulsion = 0

    # setup.xtal_length = meta.geometry.height_in_mm
    # setup.xtal_radius = meta.geometry.height_in_mm

    # setup.pc_length = meta.geometry.contact.depth_in_mm
    # setup.pc_radius = meta.geometry.contact.radius_in_mm

    # setup.wrap_around_radius = meta.geometry.groove.outer_radius_in_mm
    # setup.ditch_depth = meta.geometry.groove.depth_in_mm
    # setup.ditch_thickness = meta.geometry.groove.width_in_mm

    # setup.hole_length = setup.xtal_length - meta.geometry.well.gap_in_mm
    # # setup.hole_length_gap = meta.geometry.well.gap_in_mm
    # setup.hole_radius = meta.geometry.well.radius_in_mm

    # setup.bottom_taper_length = meta.geometry.taper.bottom.outer.height_in_mm
    # setup.outer_taper_length = meta.geometry.taper.top.outer.height_in_mm
    # setup.inner_taper_length = meta.geometry.taper.top.inner.height_in_mm
    # setup.taper_angle = meta.geometry.taper.top.inner.angle_in_deg # still confused whether this angle means outer or inner taper

    # setup.Li_thickness = meta.geometry.dl_thickness_in_mm

    # setup.xtal_HV = meta.characterization.manufacturer.op_voltage_in_V
    # setup.xtal_temp = siggen_sim.crystal_t

    # setup
# end
function detector_config(meta::PropDict, siggen_sim::SiggenSimulator)
    # function detector_config(meta::PropDict, siggen_sim::SiggenSimulator)
    siggen_geom = siggen_dict(meta, siggen_sim)
    # construct lines for siggen config file 
    detector_lines = Vector{String}()
    for (param, value) in siggen_geom
        push!(detector_lines, "$param   $value")
    end


    # read lines from user input for fieldgen 
    fieldgen_lines = readlines(open(siggen_sim.fieldgen_input))
    fieldgen_lines = replace.(fieldgen_lines, "\t" => "  ")

    total_lines = vcat(["# detector related parameters"], detector_lines, fieldgen_lines)

    # write a config file 
    sig_conf_name = siggen_sim.fieldgen_input * ".temp"
    DelimitedFiles.writedlm(sig_conf_name, total_lines)

    # construct setup
    SigGenSetup(sig_conf_name)    
end


function siggen_dict(meta::PropDict, siggen_sim::SiggenSimulator)
    Dict(
    # # hack for SSD
    #     "det_json" => meta.detector,
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
        "xtal_HV" => meta.characterization.manufacturer.op_voltage_in_V,

    # configuration for signal calculation 
        # crystal temperature in Kelvin
        "xtal_temp" => siggen_sim.crystal_t
    )
end





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

    # det_h5 = joinpath(det_path, det_name*".h5f")
    det_h5 = joinpath("cache", basename(detector_config["name"]*".h5f"))

    if !ispath(dirname(det_h5)) mkpath(dirname(det_h5)) end
    SolidStateDetectors.ssd_write(det_h5, simulation)

    @info("-> Saved cached simulation to $det_h5")

    simulation
end