abstract type PSSimulator end
struct SSDSimulator <: PSSimulator end
struct SiggenSimulator <: PSSimulator end

T = Float32

Random.seed!(123) # only for testing

# for SSD baseline and tail
# has to be larger than number of samples resulting from the online trigger filter of the DAQ
# n_baseline_samples = 2000; # SSD example was 1200;
# has to be greater than DAQ waveform length (daq_nsamples in mcpss_to_t1pss.jl)
# total_waveform_length = 8000;

function mcstp_to_mcpss(det_path::AbstractString, det_name::AbstractString, mc_events::Table, sim_config_file::AbstractString)
    sim_config = LegendGeSim.load_config(sim_config_file)
    mcstp_to_mcpss(det_path, det_name, mc_events, sim_config)

end


function mcstp_to_mcpss(det_path::AbstractString, det_name::AbstractString, mc_events::Table, sim_config::PropDict)
    noise_model = NoiseModel(sim_config)
    mcstp_to_mcpss(det_path, det_name, mc_events, noise_model)
end


##
"""
    cstp_to_mcpss(det_path, det_name, mc_events)

Simulate waveforms based on given MC events for a given detector.
Returns the resulting simulation and MC truth

det_path: path to detector json file
det_name: detector name (e.g. "V05266A")
The code will look for a .json file "det_path/det_name.json"
mc_events: table in the mcstp format (output of g4_to_mcstp)

Output: Table, Table
"""
function mcstp_to_mcpss(det_path::AbstractString, det_name::AbstractString, mc_events::Table, noise_model::NoiseModel)

    simulation = detector_simulation(det_path, det_name)

    # add fano noise, don't add if data noise is applied later
    mc_events = add_noise(mc_events, simulation, noise_model)

    # simulate waveforms
    mcpss_table, mcpss_mctruth = simulate_wf(mc_events, simulation, SSDSimulator())

    mcpss_table, mcpss_mctruth

end


"""
    detector_simulation(det_path, det_name)

Read cached h5 detector simulation if exists,
otherwise simulate detector based on json geometry.

det_path: path to detector json file
det_name: detector name (e.g. "V05266A").

Output: SSD detector simulation object
"""
function detector_simulation(det_path::AbstractString, det_name::AbstractString)

    det_h5 = joinpath(det_path, "cache", det_name*".h5f")
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
    det_h5 = joinpath(det_path, "cache", det_name*".h5f")

    if !ispath(dirname(det_h5)) mkpath(dirname(det_h5)) end
    SolidStateDetectors.ssd_write(det_h5, simulation)

    @info "Detector simulation complete"

    simulation
end


"""
    add_fnoise(mc_events, simulation, noise_model)

Add fano noise to MC events


Output: no output, the function modifies the field mc_events
"""
function add_noise(mc_events::Table, simulation::SolidStateDetectors.Simulation, ::NoiseSim)
    @info("Adding fano noise")
    det_material = simulation.detector.semiconductors[1].material
    add_fano_noise(mc_events, det_material.E_ionisation, det_material.f_fano)
end


function add_noise(mc_events::Table, ::SolidStateDetectors.Simulation, ::NoiseData)
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
