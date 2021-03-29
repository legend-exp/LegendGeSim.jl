"""
    g4_to_mcstp(g4_filename)

Process a raw g4simple simulation output (hdf5) into a Table grouped grouped by event,
which serves as input to SolidStateDetectors.simulate_waveforms() function
for waveform generation from monte-carlo simulated events.
Note: should be a part of the LegendHDF5IO package in the future.

g4_filename: path to the g4simple output hdf5 file
Output: Table
"""
function g4_to_mcstp(g4_filename::AbstractString)

    # Read raw g4simple HDF5 output file and construct arrays of corresponding data
    @info("Processing file: $g4_filename")

    ### ---- replace with LegendHDF5IO function
    g4sfile = h5open(g4_filename, "r")
    g4sntuple = g4sfile["default_ntuples"]["g4sntuple"]

    evtno = read(g4sntuple["event"]["pages"])
    detno = read(g4sntuple["iRep"]["pages"]) # no such thing for HADES
    thit = read(g4sntuple["t"]["pages"]).*u"ns"
    edep = read(g4sntuple["Edep"]["pages"]).*u"MeV"
    ekin = read(g4sntuple["KE"]["pages"]).*u"MeV"
    volID = read(g4sntuple["volID"]["pages"])
    stp = read(g4sntuple["step"]["pages"])
    mom = read(g4sntuple["parentID"]["pages"])
    trk = read(g4sntuple["trackID"]["pages"])
    pdg = read(g4sntuple["pid"]["pages"])

    x = read(g4sntuple["x"]["pages"])
    y = read(g4sntuple["y"]["pages"])
    z = read(g4sntuple["z"]["pages"])
    ### ---- replace with LegendHDF5IO function

    # Construct array of positions for input to SSD
    n_ind = length(evtno)
    pos = [ SVector{3}(([ x[i], y[i], z[i] ] .* u"mm")...) for i in 1:n_ind ]

    # println("...constructing Table")
    # Construct a Julia DataFrame with the arrays we just constructed from the g4sfile data to make grouping easier
    # raw_df = TypedTables.Table(
    raw_df = DataFrame(
            evtno = evtno,
            detno = detno,
            thit = thit,
            edep = edep,
            pos = pos,
            ekin = ekin,
            volID = volID,
            stp = stp,
            mom = mom,
            trk = trk,
            pdg = pdg
        )


    println("...group by volume")
    # Save only events that occur in the detector PV
    gdf = DataFrames.groupby(raw_df, :volID)
    # volID = 1 for the detectors in CAGE g4simple sims, this selects only events in the detector
    det_hits = DataFrame(gdf[1])

    # TODO: use Table instead of DataFrame
    # gdf = group_by_column(raw_df, :volID)
    # det_hits = gdf[1]

    # Need to turn DataFrame into a Table before using internal SSD functions (group_by_evtno, cluster_detector_hits, etc)
    # Only include parameters needed by SSD

    # println("...constructing table")
    hits_flat = Table(
        evtno = det_hits.evtno,
        detno = det_hits.detno,
        thit = det_hits.thit,
        edep = det_hits.edep,
        pos = det_hits.pos
     )

     # group hits by event number and cluster hits based on distance
    hits_by_evtno = RadiationDetectorSignals.group_by_evtno(hits_flat)
    @info("$(sum(length.(hits_by_evtno.edep))) hits before clustering")
    mc_events_clustered = @time SolidStateDetectors.cluster_detector_hits(hits_by_evtno, 0.2u"mm")
    @info("$(sum(length.(mc_events_clustered.edep))) hits after clustering")

    # Waveform generation has to be per detector.
    # Let's reshuffle the detector hits, grouping by event number and detector:
    # println("...group by detector")
    hits = RadiationDetectorSignals.ungroup_by_evtno(mc_events_clustered)
    mc_events_per_det = RadiationDetectorSignals.group_by_evtno_and_detno(hits)

    # The hits are now grouped by event number, but separately for each detector, and sorted by detector number:
    issorted(mc_events_per_det.detno)

    #This makes it easy to group them by detector number ...
    # currently there is only 1 detector with ID = 0
    mc_events_det1 = filter(evt -> evt.detno == 0, mc_events_per_det)

    mc_events_det1

end


# function group_by_column(table::TypedTables.Table, colname::Symbol)
#     sorting_idxs = sortperm(getproperty(table, colname))
#     sorted = table[sorting_idxs]
#     TypedTables.Table(consgroupedview((getproperty(table, colname), Tables.columns(sorted))))
# end