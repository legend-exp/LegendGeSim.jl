T = Float32

Random.seed!(123) # only for testing

# for SSD baseline and tail
# has to be larger than number of samples resulting from the online trigger filter of the DAQ
n_baseline_samples = 2000; # SSD example was 1200;
# has to be greater than DAQ waveform length (daq_nsamples in mcpss_to_t1pss.jl)
total_waveform_length = 8000;

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
function mcstp_to_mcpss(det_path::AbstractString, det_name::AbstractString, mc_events::Table)

    simulation = detector_simulation(det_path, det_name)

    mcstp = MCstp(simulation, mc_events)

    # add fano noise
    add_fnoise(mcstp)

    # simulate waveforms
    mcpss_table, mcpss_mctruth = simulate_wf(mcstp)

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
    read_mcstp(filename)

Helper function to read MC events from an HDF5 file with mcstp format
(g4simple events in the format compatible with SSD simulate_waveforms method)
(to be rewritten as a g4_to_mcstp() function that reads from storage)

filename: name of mcstp file
Output: Table
"""
function read_mcstp(mcfilename::AbstractString)
    # Let's read in some Monte-Carlo events (produced by Geant4).
    # We'll either read from Geant4 CSV and cache the result as HDF5,
    # or read directly from HDF5 if already available:

    mctruth_filename_csv = "data/$(mcfilename)_mcstp.csv"
    mctruth_filename_hdf5 = "cache/$(mcfilename)_mcstp.h5"

    if isfile(mctruth_filename_hdf5)
        @info("Reading MC events from HDF5.")
        mc_events = HDF5.h5open(mctruth_filename_hdf5, "r") do input
            readdata(input, "mctruth")
        end
    else
        @info("Reading MC events from Geant4-CSV.")
        mc_events = open(read, mctruth_filename_csv, Geant4CSVInput)
        mkpath(dirname(mctruth_filename_hdf5))
        HDF5.h5open(mctruth_filename_hdf5, "w") do output
            writedata(output, "mctruth", mc_events)
        end
    end

    mc_events

end

##
"""
    MCstp(simulation, mc_events)
A struct containing two things needed for waveform simulation:
    simulation: SSD detector simulation object
    mc_events: Table of mcstp format
"""
mutable struct MCstp
    simulation::SolidStateDetectors.Simulation
    mc_events::Table
end

"""
    add_fnoise(mcstp)

Add fano noise to MC events

mcstp: MCstp struct

Output: no output, the function modifies the field mc_events
"""
function add_fnoise(mcstp::MCstp)
    @info("...adding fano noise")
    det_material = mcstp.simulation.detector.semiconductors[1].material
    mcstp.mc_events = add_fano_noise(mcstp.mc_events, det_material.E_ionisation, det_material.f_fano)
end

"""
    simulate_wf(mcstp)

Simulate waveforms based on given MC events and detector simulation,
contained in the MCstp struct. Returns resulting mcpss table
and a secon table containing MC truth

mcstp: MCstp struct

Output: Table, Table
"""
function simulate_wf(mcstp::MCstp)
    # We need to filter out the few events that,
    # due to numerical effects, lie outside of the detector
    # (the proper solution is to shift them slightly, this feature will be added
    # in the future):

    filtered_events = mcstp.mc_events[findall(pts -> all(p -> T.(ustrip.(uconvert.(u"m", p))) in mcstp.simulation.detector, pts), mcstp.mc_events.pos)];

    @info("Simulating waveforms")
    contact_charge_signals = SolidStateDetectors.simulate_waveforms(
            # filtered_events[1:2000],
            # filtered_events[1:800],
            filtered_events,
            mcstp.simulation,
            max_nsteps = 4000,
            Δt = 1u"ns",
            verbose = false);

    waveforms = ArrayOfRDWaveforms(contact_charge_signals.waveform)

    # extend tail
    @info("...extending tail -> $(n_baseline_samples) baseline samples, wf length $(total_waveform_length)")
    waveforms = ArrayOfRDWaveforms(SolidStateDetectors.add_baseline_and_extend_tail.(waveforms, n_baseline_samples, total_waveform_length));

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
