Random.seed!(123) # only for testing


"""
    mstp_to_mcpss(mc_events, sim_config_file)

Simulate waveforms based on given MC events for a given configuration.
Returns the resulting simulated waveforms and MC truth

mc_events: table in the mcstp format (output of g4_to_mcstp)
sim_config_file: simulation configuration json file
"""
function mcstp_to_mcpss(mc_events::Table, sim_config_file::AbstractString)
    sim_config = LegendGeSim.load_config(sim_config_file)
    mcstp_to_mcpss(mc_events, sim_config)

end


"""
    mstp_to_mcpss(mc_events, sim_config)

Simulate waveforms based on given MC events for a given configuration.
Returns the resulting simulated waveforms and MC truth

mc_events: table in the mcstp format (output of g4_to_mcstp)
sim_config: PropDict object with simulation configuration
"""
function mcstp_to_mcpss(mc_events::Table, sim_config::PropDict)
    noise_model = NoiseModel(sim_config)
    simulation = detector_simulation(sim_config.detector_path, sim_config.detector)
    ps_simulator = PSSimulator(sim_config)
    mcstp_to_mcpss(mc_events, simulation, ps_simulator, noise_model)
end


"""
    mstp_to_mcpss(mc_events, simulation, noise_model)

Simulate waveforms based on given MC events, detector simulation and noise model.
Returns the resulting simulated waveforms and MC truth

mc_events: table in the mcstp format (output of g4_to_mcstp)
simulation: SSD simulation object
noise_model: NoiseModel object

Output: Table, Table
"""
function mcstp_to_mcpss(mc_events::Table, simulation::SolidStateDetectors.Simulation, ps_simulator::PSSimulator, noise_model::NoiseModel)
    # add fano noise, don't add if data noise is applied later
    mc_events = add_fano_noise(mc_events, simulation, noise_model)

    # simulate waveforms
    mcpss_table, mcpss_mctruth = simulate_wf(mc_events, simulation, ps_simulator)

    mcpss_table, mcpss_mctruth
end


"""
    detector_simulation(det_path, det_name)

Read cached h5 detector simulation if exists,
otherwise simulate detector based on json geometry.

det_path: path to detector json file
det_name: detector name without extension (e.g. "IC160A")

Output: SSD detector simulation object
"""
function detector_simulation(det_path::AbstractString, det_name::AbstractString)

    det_h5 = joinpath(det_path, det_name*".h5f")
    if isfile(det_h5)
        @info "Reading $det_name simulation from cached h5"
        simulation = SolidStateDetectors.ssd_read(det_h5, Simulation)
    else
        @info "Simulating $det_name from scratch"
        simulation = simulate_detector(det_path, det_name)
    end

    simulation
end


"""
    simulate_detector(det_path, det_name)

Simulate detector based on json geometry.

det_path: path to detector json file
det_name: detector name (e.g. "V05266A").
The code will look for a .json file "det_path/det_name.json"

Output: SSD detector simulation object
"""
function simulate_detector(det_path::AbstractString, det_name::AbstractString)

    det_geom = joinpath(det_path, det_name * ".json")
    @info("Reading geometry from $det_geom")

    simulation = Simulation{T}(det_geom)

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

    @info("-> Saving to cache/.h5...")
    det_h5 = joinpath(det_path, det_name*".h5f")

    if !ispath(dirname(det_h5)) mkpath(dirname(det_h5)) end
    SolidStateDetectors.ssd_write(det_h5, simulation)

    @info "Detector simulation complete"

    simulation
end


