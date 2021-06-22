"""
    pet_to_stp(pet_filename, detector_SSD)

AbstractString, SolidStateDetectors.Simulation -> Table    

Construct a table with "stepping info" based on
    position-energy-time information of simulated energy depositions
    given in <pet_filename> and detector geometry given in <detector_SSD>.

Current steps include:
    - removing events outside of detector PV
    - clustering
    - grouping by detector
    - removing events reconstructed to be outisde of the detector
"""
function pet_to_stp(pet_filename::AbstractString, detector_SSD::SolidStateDetectors.Simulation)
    @info "---------------------- pet -> stp (stepping info)"

    # Read raw g4simple HDF5 output file and construct arrays of corresponding data
    println("Processing file: $pet_filename")

    # currently PET info is simulated by g4simple
    pet_table = open(pet_filename, Geant4HDF5Input) do io
        read(io)
    end
    # pet_table = read_pet(pet_filename)

    # Save only events that occur in the detector PV
    # volID = 1 selects only events in the detector
    pet_pv = group_by_column(pet_table, :volID)[1]

    # Need to turn into a normal Table before using internal SSD functions (group_by_evtno, cluster_detector_hits, etc)
    pet_flat = Table(
        evtno = pet_pv.evtno,
        detno = pet_pv.detno,
        thit = pet_pv.thit,
        edep = pet_pv.edep,
        pos = pet_pv.pos
    )

    @info "...clustering"
    # group hits by event number and cluster hits based on distance
    hits_by_evtno = RadiationDetectorSignals.group_by_evtno(pet_flat)
    println("$(sum(length.(hits_by_evtno.edep))) hits before clustering")
    events_clustered = @time SolidStateDetectors.cluster_detector_hits(hits_by_evtno, 0.2u"mm")
    println("$(sum(length.(events_clustered.edep))) hits after clustering")

    # Waveform generation has to be per detector.
    # Let's reshuffle the detector hits, grouping by event number and detector
    @info("...grouping by detector")
    hits = RadiationDetectorSignals.ungroup_by_evtno(events_clustered)
    events_per_det = RadiationDetectorSignals.group_by_evtno_and_detno(hits)

    # The hits are now grouped by event number, but separately for each detector, and sorted by detector number:
    issorted(events_per_det.detno) # what does this do? It's not "!"

    # This makes it easy to group them by detector number ...
    # currently there is only 1 detector with ID = 0
    events_det1 = filter(evt -> evt.detno == 0, events_per_det)

    @info "...removing events outside of the detector"
    
    # We need to filter out the few events that, due to numerical effects, lie outside of the detector
    # (the proper solution is to shift them slightly, this feature will be added in the future)
    
    # construct SSD detector with arbitrary environment for quick event filering
    # det_config_SSD = ssd_config(det_meta, Environment())
    # detector_SSD = Simulation(SolidStateDetector{T}(det_config_SSD))

    stp_events = events_det1[findall(pts -> all(p -> T.(ustrip.(uconvert.(u"m", p))) in detector_SSD.detector, pts), events_det1.pos)]

    stp_events

end


function pet_to_stp(pet_filename::AbstractString, json_filepath::AbstractString)
    # load prop dict 
    pdict = propdict(json_filepath)

    # if the prop dict is LEGEND metadata...
    if haskey(pdict, :production)
        # ...construct SSD detector and proceed to pet->stp
        det_config_SSD = ssd_config(pdict, Environment())
        detector_SSD = Simulation(SolidStateDetector{T}(det_config_SSD))        
        pet_to_stp(pet_filename, detector_SSD)    
    # ... otherwise if the prop dict is simulation config... 
    elseif haskey(pdict, :detector_metadata)
        # ... recursively call this function with LEGEND metadata filepath
        pet_to_stp(pet_filename, pdict.detector_metadata)
    # if it's neither, complain        
    else
        @info "The given json file $(json_filepath) is neither a LEGEND metadata file
            nor a LegendGeSim simulation config!"
    end
end


# """
# AbstractString -> Table
# """
# function read_pet(pet_filename::AbstractString)
#     # check file extension and determine the type
#     file_ext = splitext(pet_filename)[end]
#     file_type = file_ext == ".csv" ? LegendTextIO.Geant4CSVInput : LegendHDF5IO.Geant4HDF5Input
    
#     open(pet_filename, file_type) do io
#         read(io)
#     end
# end



"""
    group_by_column(table, colname)

TypedTables.Table, Symbol -> TypedTables.Table 

Group table <table> by given column <colname>.

I think this function already exists somewhere.
Gotta ask Oliver and maybe use that package instead.
"""
function group_by_column(table::TypedTables.Table, colname::Symbol)
    sorting_idxs = sortperm(getproperty(table, colname))
    sorted = table[sorting_idxs]
    TypedTables.Table(consgroupedview(getproperty(sorted, colname), Tables.columns(sorted)))
end