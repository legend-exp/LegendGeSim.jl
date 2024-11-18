function simulate_fields(detector_metadata::PropDict, environment_settings::PropDict, simulation_settings::PropDict; overwrite::Bool = false)    
    env = Environment(environment_settings)
    simulator = PSSimulator(simulation_settings)

    simulate_detector(detector_metadata, env, simulator; overwrite)
end

# user generates detector metadata in code
function simulate_fields(detector_metadata::AbstractString, environment_settings::PropDict, simulation_settings::PropDict; overwrite::Bool = false)    
    simulate_fields(propdict(detector_metadata), environment_settings, simulation_settings; overwrite)
end

# user launches directly inputting separate dicts
function simulate_fields(detector_metadata::Union{AbstractString, PropDict}, environment_settings::Dict,
    simulation_settings::Dict; overwrite::Bool = false)    
    simulate_fields(detector_metadata, PropDict(environment_settings), PropDict(simulation_settings); overwrite)
end

# user launches with all settings in one dict
function simulate_fields(detector_metadata::Union{AbstractString, PropDict}, all_settings::Union{Dict,PropDict}; overwrite::Bool = false)
    simulate_fields(detector_metadata, PropDict(all_settings).environment, PropDict(all_settings).simulation; overwrite)
end

# user launches with all settings in json
function simulate_fields(detector_metadata::Union{AbstractString, PropDict}, all_settings::AbstractString; overwrite::Bool = false)
    simulate_fields(detector_metadata, propdict(all_settings); overwrite)
end



####################################
### SSD
####################################


function simulate_detector(det_meta::PropDict, env::Environment, simulator::SSDSimulator;
    overwrite::Bool = false)
    # append detector name to cached name
    # ToDo: should happen above this in simulate_fields or something, cause the same for SSD or siggen
    # filename: extract filename without extension from path
    full_cached_name = "$(det_meta.name)_$(simulator.cached_name)"

    h5fn = joinpath("cache", full_cached_name * "_ssd_fields.h5f")
    # simulate from scratch if file not found, or overwrite is asked
    # or if no cached name was given (i.e. don't cache) - such file is never saved (clumsy! ToDo)
    return if !isfile(h5fn) || overwrite
        @info("Simulating with SSD from scratch for given settings")
        # launch simulation
        sim = simulate_ssd_fields(det_meta, env, simulator)
        # cache it
        if !(simulator.cached_name == "")
            # ToDo: notify user at the beginning about cached name? So that they can see if it's wrong and stop?
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
            @info("-> Saved cached simulation to $h5fn")
        end
        sim
    else
        # read from previously cached
        sim = construct_ssd_simulation(det_meta, env, simulator)
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

function simulate_ssd_fields(det_meta::PropDict, env::Environment, simulator::SSDSimulator)
    sim = construct_ssd_simulation(det_meta, env, simulator)

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


"""
    construct_ssd_simulation(det_meta, env, sim_settings)

PropDict, Environment, SSDSimulator -> SSD.Simulation

Construct a `SolidStateDetectors.Simulation` based on geometry as given in LEGEND metadata `det_meta`
and on envorinmental settings specified in `env` and on simulational settings specified in `sim_settings`.
"""
function construct_ssd_simulation(det_meta::PropDict, env::Environment, simulator::SSDSimulator)
    T = Float32
    CS = SolidStateDetectors.Cylindrical
    sim = Simulation{T,CS}()
    sim.medium = SolidStateDetectors.material_properties[SolidStateDetectors.materials[env.medium]]
    # temporary quickfix to provide path to crystal jsons to get impurity
    # later will be read from legend-metadata similar to pylegendmeta
    # note: currently does nothing with crystal path
    sim.detector = LEGEND_SolidStateDetector(T, det_meta, env, simulator.crystal_metadata_path)

    

    if simulator.comp != "2D"
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


