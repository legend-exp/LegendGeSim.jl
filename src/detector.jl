function simulate_fields(config::LegendGeSimConfig, args...; kwargs...)
    Simulator = config.dict.simulation.method == "SSD" ? SSDSimulator : SiggenSimulator
    sim_settings = Simulator(config.dict)
    env = Environment(config)
    cached_name = config.dict.simulation.cached_name
    simulate_fields(config.dict.detector_metadata, env, sim_settings, cached_name, args...; kwargs...)
end

############
### FieldGen

function simulate_fields(
    det_meta::PropDict,
    env::Environment,
    sim_settings::SiggenSimulator,
    cached_name::AbstractString;
    overwrite::Bool = false
)
    simulate_detector(det_meta, env, cached_name, sim_settings; overwrite)
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
function simulate_detector(
        det_meta::PropDict,
        env::Environment,
        cached_name::AbstractString,
        ps_simulator::SiggenSimulator;
        overwrite::Bool = false
    )
    @info "Constructing Fieldgen/Siggen configuration file"
    # returns the name of the resulting siggen config file
    # and the name of (already or to be) generated weighting potential file
    siggen_config_name, fieldgen_wp_name = siggen_config(det_meta, env, ps_simulator, cached_name)
    # fieldgen_wp_name = joinpath("cache", fieldgen_wp_name)

    # if the WP file with such a name exists...
    if !isfile(fieldgen_wp_name) || overwrite
        #...otherwise call fieldgen
        @info "_|~|_|~|_|~|_ Fieldgen simulation"
        fieldgen(siggen_config_name)
        @info "_|~|_|~|_|~|_ Fieldgen simulation complete"
    else
        #...do nothing, siggen will later read the files based on the generated conifg
        @info "Reading cached fields from $fieldgen_wp_name"
    end

    # a SigGenSetup object
    SigGenSetup(siggen_config_name)
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

    if !ispath("cache") mkpath("cache") end

    # construct lines for future siggen config file 
    detector_lines = Vector{String}()
    push!(detector_lines, "# detector related parameters")
    for (param, value) in siggen_geom
        push!(detector_lines, "$param   $value")
    end
    println("...fieldgen configuration")
    # create filenames for future/existing fieldgen output 
    fieldgen_names = Dict(
        "drift_name" => cached_name * "_" * siggen_sim.drift_vel,
        # "drift_name" => joinpath("..", "data", "drift_vel_tcorr.tab"),
        "field_name" => cached_name * "_fieldgen_EV.dat",
        "wp_name" => cached_name * "_fieldgen_WP.dat"
    )
    if !isfile(joinpath("cache", fieldgen_names["drift_name"])) 
        default_drift_file = joinpath(dirname(@__DIR__), "examples", "configs", "drift_vel_tcorr.tab")
        cp(default_drift_file, joinpath("cache", fieldgen_names["drift_name"]))
    end
    if !isfile(joinpath("cache", fieldgen_names["field_name"])) 
        touch(joinpath("cache", fieldgen_names["field_name"]))
    end
    if !isfile(joinpath("cache", fieldgen_names["wp_name"])) 
        touch(joinpath("cache", fieldgen_names["wp_name"]))
    end

    # add corresponding lines to existing detector geometry lines
    push!(detector_lines, "")
    push!(detector_lines, "# fieldgen filenames")
    for (param, value) in fieldgen_names
        push!(detector_lines, "$param   $value")
    end

    # # read lines from user input for fieldgen
    fieldgen_lines = readlines(open(joinpath(dirname(@__DIR__), "examples", "configs", "fieldgen_settings.txt")))
    fieldgen_lines = replace.(fieldgen_lines, "\t" => "   ")

    # unite constructed detector geometry and fieldgen settings
    total_lines = vcat(detector_lines, [""], fieldgen_lines)

    # write a siggen/fieldgen config file 
    sig_conf_name = joinpath("cache", cached_name * "_siggen_config.txt")
    
    writedlm(sig_conf_name, total_lines)
    println("...fieldgen/siggen config written to $sig_conf_name")

    # return resulting fieldgen/siggen config name 
    # (currently kinda redundant but that's because we have no idea what's gonna happen later)
    sig_conf_name, fieldgen_names["wp_name"]
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



#######
### SSD

"""
    construct_ssd_simulation(det_meta::AbstractString, env::Environment, sim_settings::SSDSimulator)

Construct a `SolidStateDetectors.Simulation` based on geometry as given in LEGEND metadata `det_meta`
and on envorinmental settings specified in `env` and on simulational settings specified in `sim_settings`.

ToDo: Actually use `env`. 
"""
function construct_ssd_simulation(det_meta::PropDict, env::Environment, sim_settings::SSDSimulator)
    T = Float32
    CS = SolidStateDetectors.Cylindrical
    sim = Simulation{T,CS}()
    sim.medium = SolidStateDetectors.material_properties[SolidStateDetectors.materials[env.medium]]
    sim.detector = LEGEND_SolidStateDetector(T, det_meta)
    if sim_settings.comp != "2D"
        error("Only 2D is supported up to now.")
    end
    sim.world = begin
        crystal_radius = to_SSD_units(T, det_meta.geometry.radius_in_mm, u"mm")
        crystal_height = to_SSD_units(T, det_meta.geometry.height_in_mm, u"mm")
        SolidStateDetectors.World{T,3,CS}((
            SolidStateDetectors.SSDInterval{T,:closed,:closed,:r0,:infinite}(zero(T), crystal_radius * 1.2),
            SolidStateDetectors.SSDInterval{T,:closed,:closed,:reflecting,:reflecting}(zero(T), zero(T)),
            SolidStateDetectors.SSDInterval{T,:closed,:closed,:infinite,:infinite}(-0.2 * crystal_height, 1.2 * crystal_height)
        ))
    end
    sim.weighting_potentials = Missing[missing for i = 1:2]
    sim
end

function simulate_fields(det_meta::PropDict, env::Environment, sim_settings::SSDSimulator)
    sim = construct_ssd_simulation(det_meta, env, sim_settings)

    field_sim_settings = (
        refinement_limits = [0.2, 0.1, 0.05, 0.02, 0.01],
        depletion_handling = true,
        convergence_limit = 1e-7,
        max_n_iterations = 50000,
        verbose = false
    )

    println("...electric potential")
    calculate_electric_potential!(sim; field_sim_settings...)

    println("...electric field")
    calculate_electric_field!(sim, n_points_in_φ = 72)

    for contact in sim.detector.contacts
        println("...weighting potential $(contact.id)")
        calculate_weighting_potential!(sim, contact.id; field_sim_settings...)
    end

    return sim
end

function simulate_fields(
    det_meta::PropDict,
    env::Environment,
    sim_settings::SSDSimulator,
    cached_name::AbstractString;
    overwrite::Bool = false
)

    h5fn = joinpath("cache", cached_name * "_fields_ssd.h5f")
    return if !isfile(h5fn) || overwrite
        @info("Simulating $h5fn with SSD from scratch for given settings")
        sim = simulate_fields(det_meta, env, sim_settings)
        if !ispath(dirname(h5fn))
            mkpath(dirname(h5fn))
        end
        HDF5.h5open(h5fn, "w") do h5f
            LegendHDF5IO.writedata(h5f, "SSD_electric_potential", NamedTuple(sim.electric_potential))
            LegendHDF5IO.writedata(h5f, "SSD_point_types", NamedTuple(sim.point_types))
            LegendHDF5IO.writedata(h5f, "SSD_q_eff_fix", NamedTuple(sim.q_eff_fix))
            LegendHDF5IO.writedata(h5f, "SSD_q_eff_imp", NamedTuple(sim.q_eff_imp))
            LegendHDF5IO.writedata(h5f, "SSD_dielectric_distribution", NamedTuple(sim.ϵ_r))
            LegendHDF5IO.writedata(h5f, "SSD_electric_field", NamedTuple(sim.electric_field))
            for i in eachindex(sim.weighting_potentials)
                LegendHDF5IO.writedata(h5f, "SSD_weighting_potential_$(i)", NamedTuple(sim.weighting_potentials[i]))
            end
        end
        # SolidStateDetectors.ssd_write(h5fn, sim)
        @info("-> Saved cached simulation to $h5fn")
        sim
    else
        sim = construct_ssd_simulation(det_meta, env, sim_settings)
        @info("Reading SSD simulation from cached file $h5fn")
        HDF5.h5open(h5fn, "r") do h5f
            sim.electric_potential = ElectricPotential(LegendHDF5IO.readdata(h5f, "SSD_electric_potential"))
            sim.point_types = PointTypes(LegendHDF5IO.readdata(h5f, "SSD_point_types"))
            sim.q_eff_fix = EffectiveChargeDensity(LegendHDF5IO.readdata(h5f, "SSD_q_eff_fix"))
            sim.q_eff_imp = EffectiveChargeDensity(LegendHDF5IO.readdata(h5f, "SSD_q_eff_imp"))
            sim.ϵ_r = DielectricDistribution(LegendHDF5IO.readdata(h5f, "SSD_dielectric_distribution"))
            sim.electric_field = ElectricField(LegendHDF5IO.readdata(h5f, "SSD_electric_field"))
            for i in eachindex(sim.weighting_potentials)
                sim.weighting_potentials[i] = WeightingPotential(LegendHDF5IO.readdata(h5f, "SSD_weighting_potential_$(i)"))
            end
        end
        sim
    end
end
