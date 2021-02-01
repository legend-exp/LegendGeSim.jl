##############################################
# Written by Gulden Othman for CAGE simulation processing. Should work for any g4simple output file in hdf5 format
# July, 2020
##############################################
# https://github.com/legend-exp/CAGE/blob/master/sims/post_process/julia_post_process.jl

@info "Loading packages..."
using LegendDataTypes, LegendHDF5IO
using ArraysOfArrays, StaticArrays, Tables, TypedTables
using HDF5
using DataFrames
using SolidStateDetectors: cluster_detector_hits
using RadiationDetectorSignals: group_by_evtno, ungroup_by_evtno, group_by_evtno_and_detno
using Unitful


function main()
    g4_dir = "/lfs/l1/legend/detector_char/enr/hades/simulations/legend-g4simple-simulation/IC-legend/IC160A/Th228/uncollimated/top_source_holder/hdf5/"
#    g4_dir = "data/"
    processed_dir = "cache/"
    base_filename = "raw-IC160A-Th228-uncollimated-top-run0002-source_holder-bi-hdf5-01-test"
    raw_extension = ".hdf5"
    processed_extension = ".h5"

    @info "Processing g4 events"
    mcstp = g4_to_mcstp(g4_dir, processed_dir, base_filename, raw_extension, processed_extension, save=true)
end


function g4_to_mcstp(g4_dir, processed_dir, base_filename, raw_extension, processed_extension; save::Bool=false)
    # Use when you want to process a raw g4simple simulation output (hdf5) into a Table grouped by event
    # This produces a Table in a .lh5 file that can then be directly input into the SSD
    # SolidStateDetectors.simulate_waveforms() function for waveform generation from monte-carlo simulated events


    # Read raw g4simple HDF5 output file and construct arrays of corresponding data
    filename = g4_dir * base_filename * raw_extension
    println("Processing file: $filename")
    g4sfile = h5open(filename, "r")
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

    # Translate z-coordinates to match SSD geometry-- CAGE specific
    # if icpc==true
    #     println("Translating z-dimension so surface is at ICPC height")
    #     z .+= (22.5 + 86.4)
    # elseif oppi==true
    #     println("Translating z-dimension so surface is at OPPI height")
    #     z .+= (22. + 51.)
    # else
    #     println("Keeping in g4simple coordinates")
    # end

    # Construct array of positions for input to SSD
    n_ind = length(evtno)
    pos = [ SVector{3}(([ x[i], y[i], z[i] ] .* u"mm")...) for i in 1:n_ind ]


    println("...constructing DataFrame")
    # Construct a Julia DataFrame with the arrays we just constructed from the g4sfile data to make grouping easier
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
    gdf = groupby(raw_df, :volID)

    # volID = 1 for the detectors in CAGE g4simple sims, this selects only events in the detector
#    det_hits = DataFrame(gdf[2])
    det_hits = DataFrame(gdf[1]) # for HADES volID = 1

    # Need to turn DataFrame into a Table before using internal SSD functions (group_by_evtno, cluster_detector_hits, etc)
    # Only include parameters needed by SSD

    println("...constructing table")
    hits_flat = Table(
        evtno = det_hits.evtno,
        detno = det_hits.detno,
        thit = det_hits.thit,
        edep = det_hits.edep,
        pos = det_hits.pos
     )

     # group hits by event number and cluster hits based on distance
    hits_by_evtno = group_by_evtno(hits_flat)
#    hits_clustered = cluster_detector_hits(hits_by_evtno, 0.2u"mm")
#    hits_by_det = group_by_evtno_and_detno(ungroup_by_evtno(hits_clustered))

    # create output filename and save Table to .lh5
    out_filename = processed_dir *  base_filename * "_mcstp" * processed_extension

    # println(typeof(out_filename))

    # Save output to .lh5
    if(save)
        h5open(out_filename, "w") do f
            LegendDataTypes.writedata(f, "mctruth", hits_by_evtno)
        end
        println("Processed file save to: $out_filename")
    end

    hits_by_evtno


end

##

#main()
