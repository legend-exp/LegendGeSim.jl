"""
Decide which function to dispatch depending on whether
a LEGEND metadata json is given 
or a siggen config file (that contains detector geometry)
"""
function g4_to_mcstp(g4_filename::AbstractString, sim_config::PropDict)
    g4_to_mcstp(g4_filename, sim_config.detector)
end


function g4_to_mcstp(g4_filename::AbstractString, detector::AbstractString)
    det_config_SSD = detector_config(detector, Environment(), SSDSimulator())
    g4_to_mcstp(g4_filename, det_config_SSD)
end


"""
    g4_to_mcstp(g4_filename)

Process a raw g4simple simulation output (hdf5) into a Table grouped grouped by event,
which serves as input to SolidStateDetectors.simulate_waveforms() function
for waveform generation from monte-carlo simulated events.
Note: should be a part of the LegendHDF5IO package in the future.

g4_filename: path to the g4simple output hdf5 file
Output: Table
"""
function g4_to_mcstp(g4_filename::AbstractString, det_config_SSD::Dict)
    # Read raw g4simple HDF5 output file and construct arrays of corresponding data
    @info("Processing file: $g4_filename")

    g4_table = open(g4_filename, Geant4HDF5Input) do io
        read(io)
    end

    # Save only events that occur in the detector PV
    # volID = 1 selects only events in the detector
    g4_pv = group_by_column(g4_table, :volID)[1]

    # Need to turn into a normal Table before using internal SSD functions (group_by_evtno, cluster_detector_hits, etc)

    g4_flat = Table(
        evtno = g4_pv.evtno,
        detno = g4_pv.detno,
        thit = g4_pv.thit,
        edep = g4_pv.edep,
        pos = g4_pv.pos
    )

    # group hits by event number and cluster hits based on distance
    hits_by_evtno = RadiationDetectorSignals.group_by_evtno(g4_flat)
    @info("$(sum(length.(hits_by_evtno.edep))) hits before clustering")
    mc_events_clustered = @time SolidStateDetectors.cluster_detector_hits(hits_by_evtno, 0.2u"mm")
    @info("$(sum(length.(mc_events_clustered.edep))) hits after clustering")

    # Waveform generation has to be per detector.
    # Let's reshuffle the detector hits, grouping by event number and detector
    @info("Grouping by detector...")
    hits = RadiationDetectorSignals.ungroup_by_evtno(mc_events_clustered)
    mc_events_per_det = RadiationDetectorSignals.group_by_evtno_and_detno(hits)

    # The hits are now grouped by event number, but separately for each detector, and sorted by detector number:
    issorted(mc_events_per_det.detno) # what does this do? It's not "!"

    #This makes it easy to group them by detector number ...
    # currently there is only 1 detector with ID = 0
    mc_events_det1 = filter(evt -> evt.detno == 0, mc_events_per_det)

    # Remove events lying outside of the detector
    @info "Remove events outside of the detector..."
    detector_SSD = Simulation(SolidStateDetector{T}(det_config_SSD))
    # We need to filter out the few events that,
    # due to numerical effects, lie outside of the detector
    # (the proper solution is to shift them slightly, this feature will be added
    # in the future):
    mc_events = mc_events_det1[findall(pts -> all(p -> T.(ustrip.(uconvert.(u"m", p))) in detector_SSD.detector, pts), mc_events_det1.pos)]

    mc_events

end




function group_by_column(table::TypedTables.Table, colname::Symbol)
    sorting_idxs = sortperm(getproperty(table, colname))
    sorted = table[sorting_idxs]
    TypedTables.Table(consgroupedview(getproperty(sorted, colname), Tables.columns(sorted)))
end