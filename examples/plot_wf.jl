using LegendGeSim
using Plots
using HDF5

using LegendHDF5IO: readdata, writedata
using TypedTables

##
mc_name = "raw-IC160A-Th228-uncollimated-top-run0002-source_holder-bi-hdf5-01-test"
mc_path = "cache/"

filename = joinpath(mc_path, mc_name*"_mcpss.h5")
mcpss = LegendGeSim.read_mcpss(filename)

##

plot(mcpss.waveform[1:10])

##

function read_raw(filename, path)
    println("File: $filename")
    HDF5.h5open(filename) do input
        Table(
            baseline = readdata(input, "$path/baseline"),
            channel = readdata(input, "$path/channel"),
            energy = readdata(input, "$path/energy"),
            ievt = readdata(input, "$path/ievt"),
            numtraces = readdata(input, "$path/numtraces"),
            packet_id = readdata(input, "$path/packet_id"),
            timestamp = readdata(input, "$path/timestamp"),
            waveform = readdata(input, "$path/waveform"),
            wf_max = readdata(input, "$path/wf_max"),
            wf_std = readdata(input, "$path/wf_std"),
        )
    end
end

##

filename = joinpath(mc_path, mc_name*"_mcraw.h5")
mcraw = read_raw(filename, "raw")

##

plot(mcraw.waveform[1:10])