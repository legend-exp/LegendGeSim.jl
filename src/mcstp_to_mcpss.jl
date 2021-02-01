# Convert mcraw HDF5 format to mcpss
# 19.01.2021 Mariia Redchuk mariia.redchuk@pd.infn.it


@info "Loading packages..."

using ArraysOfArrays, StaticArrays, Tables, TypedTables
using Statistics, Random, Distributions, StatsBase
using Unitful
using RadiationDetectorSignals
using RadiationDetectorSignals: group_by_evtno, ungroup_by_evtno, group_by_evtno_and_detno
using HDF5
using SolidStateDetectors

using LegendDataTypes
using LegendHDF5IO
using LegendHDF5IO: readdata, writedata
using LegendTextIO # Geant4CSVInput

import HDF5

T = Float32

Random.seed!(123) # only for testing

# for SSD baseline and tail
# has to be larger than number of samples resulting from the online trigger filter of the DAQ
n_baseline_samples = 2000; # SSD example was 1200;
# has to be greater than DAQ waveform length (daq_nsamples in mcpss_to_t1pss.jl)
total_waveform_length = 8000;

##

function main()
    det_name = "V05266A"
#    mc_name = "raw-IC160A-Th228-uncollimated-top-run0002-source_holder-bi-hdf5-02"
    mc_name = "raw-IC160A-Th228-uncollimated-top-run0002-source_holder-bi-hdf5-01-test"

    # prepare clustered events for single detector (already in h5 or CSV)
    mc_events = read_mcstp(mc_name)
    mc_events = prepare_mcstp(mc_events)
    mcstp_to_mcpss(det_name, mc_events, mc_name=mc_name)
end

function mcstp_to_mcpss(det_name, mc_events; mc_name::String="")
    println("Detector simulation")
    if isfile("cache/$(det_name).h5f")
        simulation = SolidStateDetectors.ssd_read("cache/$(det_name).h5f", Simulation)
    else
        simulation = simulate_detector(det_name)
    end


    mcstp = MCstp(simulation, mc_events)

    # add fano noise
    add_fnoise(mcstp)

    # simulate waveforms
    mcpss_table, mcpss_mctruth = simulate_wf(mcstp)

    # save results
    if(mc_name != "")
        println("...saving table")
        out_filename = "cache/$(mc_name)_mcpss.h5"
        println("-> $out_filename")
        h5open(out_filename, "w") do f
            LegendDataTypes.writedata(f, "mcpss/mcpss", mcpss_table)
            LegendDataTypes.writedata(f, "mcpss/mctruth", mcpss_mctruth)
        end
        println("Done")
    end

    mcpss_table, mcpss_mctruth 

end

function simulate_detector(det_name)
    det_geom = "data/$(det_name).json"
    println("Reading geometry from $det_geom")

    simulation = Simulation{T}(det_geom)

    println("-> Electric potential...")
    calculate_electric_potential!( simulation,
                               max_refinements = 4)

    println("-> Electric field...")
    calculate_electric_field!(simulation, n_points_in_φ = 72)

    println("-> Capacitance...")
    calculate_capacitance(simulation)

    println("-> Drift field...")
    calculate_drift_fields!(simulation)

    println("-> Weighting potential...")
    for contact in simulation.detector.contacts
        calculate_weighting_potential!(simulation, contact.id, max_refinements = 4, n_points_in_φ = 2, verbose = false)
    end

    println("-> Saving...")
    det_h5 = "cache/$(det_name).h5f"

    if !ispath(dirname(det_h5)) mkpath(dirname(det_h5)) end
    SolidStateDetectors.ssd_write(det_h5, simulation)

    @info "Detector simulation complete"

    simulation
end


function read_mcstp(mcfilename)
    # Let's read in some Monte-Carlo events (produced by Geant4).
    # We'll either read from Geant4 CSV and cache the result as HDF5,
    # or read directly from HDF5 if already available:

    mctruth_filename_csv = "data/$(mcfilename)_mcstp.csv"
    mctruth_filename_hdf5 = "cache/$(mcfilename)_mcstp.h5"

    if isfile(mctruth_filename_hdf5)
        println("Reading MC events from HDF5.")
        mc_events = HDF5.h5open(mctruth_filename_hdf5, "r") do input
            readdata(input, "mctruth")
        end
    else
        println("Reading MC events from Geant4-CSV.")
        mc_events = open(read, mctruth_filename_csv, Geant4CSVInput)
        mkpath(dirname(mctruth_filename_hdf5))
        println("Writing MC events to HDF5.")
        HDF5.h5open(mctruth_filename_hdf5, "w") do output
            writedata(output, "mctruth", mc_events)
        end
    end

    mc_events

end

##
function prepare_mcstp(mc_events)

    # Producing pulse shapes from raw MC events is wastful,
    # it's more efficient to cluster detectors hits (within a small radius) first:
    println("$(sum(length.(mc_events.edep))) hits before clustering")
    mc_events_clustered = @time SolidStateDetectors.cluster_detector_hits(mc_events, 0.2u"mm")
    println("$(sum(length.(mc_events_clustered.edep))) hits after clustering")

    # Waveform generation has to be per detector.
    # Let's reshuffle the detector hits, grouping by event number and detector:
    println("...group by detector")
    hits = ungroup_by_evtno(mc_events_clustered)
    mc_events_per_det = group_by_evtno_and_detno(hits)

    # The hits are now grouped by event number, but separately for each detector, and sorted by detector number:
    issorted(mc_events_per_det.detno)

    #This makes it easy to group them by detector number ...
    # currently there is only 1 detector with ID = 0
    mc_events_det1 = filter(evt -> evt.detno == 0, mc_events_per_det)

    mc_events_det1


end

##

mutable struct MCstp
    simulation
    mc_events
end


function add_fnoise(mcstp)
    println("...adding fano noise")
    det_material = mcstp.simulation.detector.semiconductors[1].material
    mcstp.mc_events = add_fano_noise(mcstp.mc_events, det_material.E_ionisation, det_material.f_fano)
end


function simulate_wf(mcstp)
    # Also, we need to filter out the few events that,
    # due to numerical effects, lie outside of the detector
    # (the proper solution is to shift them slightly, this feature will be added in the future):

    filtered_events = mcstp.mc_events[findall(pts -> all(p -> T.(ustrip.(uconvert.(u"m", p))) in mcstp.simulation.detector, pts), mcstp.mc_events.pos)];

    println("Simulating waveforms")
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
    println("...extending tail -> $(n_baseline_samples) baseline samples, wf length $(total_waveform_length)")
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

##


#main()
