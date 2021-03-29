"""
Abstract type for hierarchy and multiple dispatch
that allows to choose between SolidStateDetectors and Siggen
as waveform simulation method
"""
abstract type PSSimulator end
struct SSDSimulator <: PSSimulator end
struct SiggenSimulator <: PSSimulator end

function PSSimulator(sim_config::PropDict)
    @info "Waveform simulation method: $(sim_config.sim_method)"
    if(sim_config.sim_method == "SSD")
        SSDSimulator()
    elseif(sim_config.sim_method == "siggen")
        SiggenSimulator()
    else
        println("This simulation method is not implemented!")
    end
end

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


"""
    add_fnoise(mc_events, simulation, noise_model)

Add fano noise to MC events

mc_events: table in the mcstp format (output of g4_to_mcstp)
noise_model: NoiseSim object

Output: Table
"""
function add_fano_noise(mc_events::Table, simulation::SolidStateDetectors.Simulation, ::NoiseSim)
    @info("Adding fano noise")
    det_material = simulation.detector.semiconductors[1].material
    add_fano_noise(mc_events, det_material.E_ionisation, det_material.f_fano)
end


function add_fano_noise(mc_events::Table, ::SolidStateDetectors.Simulation, ::NoiseData)
    # do nothing since if we're using noise from data we do not simulate fano noise
    # not to double count
    @info("Not adding fano noise because using noise levels from data")
    mc_events
end


"""
    simulate_wf(mc_events, simulation, simulator)

Simulate waveforms based on given MC events and detector simulation,
contained in the MCstp struct. Returns resulting mcpss table
and a secon table containing MC truth

mcstp: MCstp struct
simulator: Simulator object inheriting from PSSimulator for multiple dispatch

Output: Table, Table
"""
function simulate_wf(mc_events::Table, simulation::SolidStateDetectors.Simulation, ::SSDSimulator)
    # We need to filter out the few events that,
    # due to numerical effects, lie outside of the detector
    # (the proper solution is to shift them slightly, this feature will be added
    # in the future):

    filtered_events = mc_events[findall(pts -> all(p -> T.(ustrip.(uconvert.(u"m", p))) in simulation.detector, pts), mc_events.pos)];

    @info("Simulating waveforms")
    contact_charge_signals = SolidStateDetectors.simulate_waveforms(
            # filtered_events[1:2000],
            # filtered_events[1:800],
            filtered_events,
            simulation,
            max_nsteps = 4000,
            Δt = 1u"ns",
            verbose = false);

    waveforms = ArrayOfRDWaveforms(contact_charge_signals.waveform)

    # extend tail
    # @info("...extending tail -> $(n_baseline_samples) baseline samples, wf length $(total_waveform_length)")
    # waveforms = ArrayOfRDWaveforms(SolidStateDetectors.add_baseline_and_extend_tail.(waveforms, n_baseline_samples, total_waveform_length));

    # convert to Tier1 format
    mcpss_table = Table(
        channel = contact_charge_signals.chnid,
        ievt = contact_charge_signals.evtno,
        waveform = waveforms
    )

    mcpss_mctruth = Table(
        detno = contact_charge_signals.detno,
        edep = contact_charge_signals.edep,
        ievt = contact_charge_signals.evtno,
        pos = contact_charge_signals.pos,
        thit = contact_charge_signals.thit
    )

    mcpss_table, mcpss_mctruth
end


function simulate_wf(mc_events::Table, simulation::SolidStateDetectors.Simulation, ::SiggenSimulator)
    # to be implemented with Julia siggen wrapper
    # which method to use, SSD or siggen, should be in the configuration
    println("This method is not yet implemented")
end
