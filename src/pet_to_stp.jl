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
function pet_to_stp(pet_table::Table, detector_SSD::SolidStateDetectors.Simulation)
    ## in case it's the Geant4 CSV file from LegendTestData, it's already grouped correctly
    hits_by_evtno = pet_table

    ## in case it's a g4simple HDF5 file

    # Save only events that occur in the detector PV
    # volID = 1 selects only events in the detector
    if hasproperty(pet_table, :volID)
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

    # Waveform generation has to be per detector.
    # Let's reshuffle the detector hits, grouping by event number and detector
    @info("...grouping by detector")
    hits = RadiationDetectorSignals.ungroup_by_evtno(events_clustered)
    events_per_det = RadiationDetectorSignals.group_by_evtno_and_detno(hits)

    # The hits are now grouped by event number, but separately for each detector, and sorted by detector number:
    issorted(events_per_det.detno) # what does this do? It's not "!"

    # This makes it easy to group them by detector number ...
    # currently there is only 1 detector with ID = 0; in Geant4 CSV in TestData it's 1
    events_det1 = filter(evt -> evt.detno == 1, events_per_det)

    @info "...removing events outside of the detector"
    
    # We need to filter out the few events that, due to numerical effects, lie outside of the detector
    # (the proper solution is to shift them slightly, this feature will be added in the future)
    stp_events = events_det1[findall(pts -> all(p -> T.(ustrip.(uconvert.(u"m", p))) in detector_SSD.detector, pts), events_det1.pos)]

    stp_events

end



function pet_to_stp(pet_filename::AbstractString, det_metadata::AbstractString)
    # construct simplified config based on given info
    sim_config = PropDict(:detector_metadata => det_metadata, :input_file => pet_filename)

    pet_to_stp(sim_config)
end


function pet_to_stp(sim_config_file::AbstractString)
    # given config file already contains inputs as well
    pet_to_stp(propdict(sim_config_file))
end


function pet_to_stp(sim_config::PropDict)
    @info "---------------------- pet -> stp (stepping info)"
    # SSD detector for geometry purposes (not simulation)
    det_config_SSD = ssd_config(sim_config.detector_metadata, Environment())
    detector_SSD = Simulation(SolidStateDetector{T}(det_config_SSD))        

    println("Processing file: $(sim_config.input_file)")    
    pet_table = read_pet(sim_config.input_file)

    # launch pet->stp 
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


# function read_pet(pet_filename::AbstractString, file_type::Type{Geant4CSVInput})
#     open(pet_filename, file_type) do io
#         read(io)
#     end
# end    


# function read_pet(pet_filename::AbstractString, ::Type{Geant4HDF5Input})
#     g4sfile = h5open(pet_filename, "r")
#     g4sntuple = g4sfile["default_ntuples"]["g4sntuple"]

#     evtno = read(g4sntuple["event"]["pages"])
#     detno = read(g4sntuple["iRep"]["pages"])
#     thit = read(g4sntuple["t"]["pages"]).*u"ns"
#     edep = read(g4sntuple["Edep"]["pages"]).*u"MeV"
#     volID = read(g4sntuple["volID"]["pages"])

#     x = read(g4sntuple["x"]["pages"])
#     y = read(g4sntuple["y"]["pages"])
#     z = read(g4sntuple["z"]["pages"])

#     # Construct array of positions for input to SSD
#     n_ind = length(evtno)
#     pos = [ SVector{3}(([ x[i], y[i], z[i] ] .* u"mm")...) for i in 1:n_ind ]

#     # Construct a Table with the arrays we just constructed from the g4sfile data
#     TypedTables.Table(
#             evtno = evtno,
#             detno = detno,
#             thit = thit,
#             edep = edep,
#             pos = pos,
#             volID = volID,
#     )
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
