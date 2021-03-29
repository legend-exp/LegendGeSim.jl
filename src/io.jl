"""
    read_mcstp(filename)

Helper function to read MC events from an HDF5 file with mcstp format
(g4simple events after clustering in the format compatible with SSD simulate_waveforms method)

To be rewritten as a g4_to_mcstp() function that reads from storage?
To be moved to LegendHDF5IO?

filename: full path to mcstp file
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



"""
    read_mcpss(filename)

Helper function to read simulation output from an HDF5 file with mcpss format
(simulated waveforms)


To be rewritten as a g4_to_mcstp() function that reads from storage?
To be moved to LegendHDF5IO?

filename: name of mcpss file
Output: Table
"""
function read_mcpss(filename::AbstractString)
    @info "Reading mcpss from $filename"

    HDF5.h5open(filename) do input
        Table(
            channel = LegendDataTypes.readdata(input, "mcpss/mcpss/channel"),
            ievt = LegendDataTypes.readdata(input, "mcpss/mcpss/ievt"),
            waveform = LegendDataTypes.readdata(input, "mcpss/mcpss/waveform")
        )
    end
end


"""
    read_mctruth(filename)

Helper function to read mctruth from an HDF5 file with mcpss format
(contains mc truth in a dedicated group)

To be rewritten as a g4_to_mcstp() function that reads from storage?
To be moved to LegendHDF5IO?

filename: name of mcpss file
Output: Table
"""
function read_mctruth(filename::AbstractString)
    @info "Reading MC truth from $filename"

    HDF5.h5open(filename) do input
        Table(
            detno = LegendDataTypes.readdata(input, "mcpss/mctruth/detno"),
            edep = LegendDataTypes.readdata(input, "mcpss/mctruth/edep"),
            ievt = LegendDataTypes.readdata(input, "mcpss/mctruth/ievt"),
            pos = LegendDataTypes.readdata(input, "mcpss/mctruth/pos"),
            thit = LegendDataTypes.readdata(input, "mcpss/mctruth/thit")
        )
    end
end