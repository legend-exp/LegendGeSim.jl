"""
    pet_to_stp(pet_table, detector_SSD)

    Table, SolidStateDetectors.Simulation -> Table    

Construct a table with "stepping info" based on
    position-energy-time information of simulated energy depositions
    given in <pet_table> and detector geometry given in <detector_SSD>.

Current steps include:
    - removing events outside of detector PV
    - clustering
    - grouping by detector
    - removing events reconstructed to be outisde of the detector
"""
function pet_to_stp(pet_table::Table, detector_SSD::SolidStateDetector)
    ## in case it's the Geant4 CSV file from LegendTestData, it's already grouped correctly
    hits_by_evtno = pet_table

    ## in case it's a g4simple HDF5 file

    # Save only events that occur in the detector PV
    # volID = 1 selects only events in the detector
    if hasproperty(pet_table, :volID)
        @info "...checking volID"
        pet_pv = group_by_column(pet_table, :volID)[1]

        # Need to turn into a normal Table before using internal SSD functions (group_by_evtno, cluster_detector_hits, etc)
        pet_flat = Table(
            evtno = pet_pv.evtno,
            detno = pet_pv.detno,
            thit = pet_pv.thit,
            edep = pet_pv.edep,
            pos = pet_pv.pos
        )

        # group hits by event number 
        hits_by_evtno = RadiationDetectorSignals.group_by_evtno(pet_flat)
    end

    # cluster hits based on distance
    @info "...clustering"
    println("$(sum(length.(hits_by_evtno.edep))) hits before clustering")
    events_clustered = @time SolidStateDetectors.cluster_detector_hits(hits_by_evtno, 0.2u"mm")
    println("$(sum(length.(events_clustered.edep))) hits after clustering")

    hits = RadiationDetectorSignals.ungroup_by_evtno(events_clustered)

    # Waveform generation has to be per detector.
    # Let's reshuffle the detector hits, grouping by event number and detector
    # @info("...grouping by detector")
    events_per_det = RadiationDetectorSignals.group_by_evtno_and_detno(hits)

    # The hits are now grouped by event number, but separately for each detector, and sorted by detector number:
    # issorted(events_per_det.detno) # what does this do? It's not "!"    

    ## ---------- This part is done for the case of multiple detectors
    ## For now we consider that the input file should have only one detector
    ## TODO
    # This makes it easy to group them by detector number ...
    # events_det1 = filter(evt -> evt.detno == 0, events_per_det)
    ## ---------- This part is done for the case of multiple detectors

    @info "...removing events outside of the detector"

    # We need to filter out the few events that, due to numerical effects, lie outside of the detector
    # (the proper solution is to shift them slightly, this feature will be added in the future)
    T = Float32
    stp_events = events_per_det[findall(pts -> all(p -> CartesianPoint{T}(T.(ustrip.(uconvert.(u"m", p)))) in detector_SSD, pts), events_per_det.pos)]
    # stp_events = events_det1[findall(pts -> all(p -> T.(ustrip.(uconvert.(u"m", p))) in detector_SSD.detector, pts), events_det1.pos)]
    # stp_events = if S isa SolidStateDetectors.Simulation
    #     events_per_det[findall(pts -> all(p -> CartesianPoint{T}(T.(ustrip.(uconvert.(u"m", p)))) in detector_SSD.detector, pts), events_per_det.pos)]
    # else
    #     events_per_det
    # end
    stp_events
end



function pet_to_stp(pet_filename::AbstractString, sim_config::LegendGeSimConfig)

    @info "---------------------- pet -> stp (stepping info)"

    # SSD detector for geometry purposes (not simulation)
    detector_SSD = LEGEND_SolidStateDetector(Float32, sim_config.dict.detector_metadata)
    pet_table = read_pet(pet_filename)

    println("Processing file: $(pet_filename) for detector $(sim_config.dict.detector_metadata.det_name)")    
    
    pet_to_stp(pet_table, detector_SSD)
end


"""
AbstractString -> Table
"""
function read_pet(pet_filename::AbstractString)
    # check file extension and determine the type
    file_ext = splitext(pet_filename)[end]
    file_type = file_ext == ".csv" ? LegendTextIO.Geant4CSVInput : LegendHDF5IO.Geant4HDF5Input
    
    # read_pet(pet_filename, file_type)
    open(pet_filename, file_type) do io
        read(io)
    end
end


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
