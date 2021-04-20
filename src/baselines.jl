# ##
# tier1_path = "/lfs/l1/legend/detector_char/enr/hades/char_data/V05266A/tier1/th_HS2_lat_psa/pygama/v01.00"
# tier1_file = "char_data-V05266A-th_HS2_lat_psa-run0001-200904T143805_tier1.lh5"

# ##

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
            # here, readdata will read the waveform in the format of RDWaveform
            waveform = readdata(input, "$path/waveform"),
            wf_max = readdata(input, "$path/wf_max"),
            wf_std = readdata(input, "$path/wf_std"),
        )
    end
end

##

# # after we read this function, we obtain a table
# # one of its columns is "waveform" and it contains RDWaveform objects
# filename = joinpath(tier1_path, tier1_file)
# tier1_table = read_raw(filename, "raw")
# ##
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
function risepoint(wf::RDWaveform)
    #define the index of wf rising point
    base_start = wf.value[1] 
    risepoint_index = 0
    for i in 1:1950
        if wf.value[i] < base_start + 5
            risepoint_index = i
        end
    end
    return risepoint_index
end

function basestart(wfs::ArrayOfRDWaveforms)       
    #get the distribution of baseline starting point 
    basestart_list = []
    for i in wfs
        append!(basestart_list, i.value[1])
    end
    uplim = mean(basestart_list) + 1000
    dolim = mean(basestart_list) - 1000
    return uplim, dolim
end  

function peak_finder(wf::RDWaveform,n)
    peak = findmax(wf.value)[2]
    return RDWaveform(wf.time[peak:peak+n], wf.value[peak:peak+n])
end

function extract_baseline(wf::RDWaveform)
    #extract baseline
    baseline_index = risepoint(wf)
    return RDWaveform(wf.time[begin:baseline_index], wf.value[begin:baseline_index])
end

function base_slope(wf::RDWaveform, n)
    #slope of the tail  n = 1800
peak = findmax(wf.value)[2]
charge = wf.value[begin:peak-n]
time = wf.time[begin:peak-n]
f = linear_fit(ustrip(time), ustrip(charge))
slope = f[2]
return slope
end

function tail_slope(wf::RDWaveform, n)
    #slope of the tail, n = 500, n= 1500, n=2000
peak = findmax(wf.value)[2]
charge = wf.value[peak:peak+n]
time = wf.time[peak:peak+n]
f = linear_fit(ustrip(time), ustrip(charge))
slope = f[2]
return slope
end


function stats(sample, accepted_list, tau_list, stdlist)
    #stats
println("No. of waveforms = ", sample) 
println("Accepted = ", length(accepted_list), "  ----------  ", (length(accepted_list)/sample)*100 , " %") 
println("mean std = ", mean(stdlist))                                  #means std
println("max tau = ", -1/maximum(tau_list)*16/1000 , " microsecond")    #max tau
println("min tau = ", -1/minimum(tau_list)*16/1000 , " microsecond")    #min tau
end

function selection_cut(wf::RDWaveform, sample)
    #main
    std_list = []               # store std for every wf
    tau_list = []               # store slope of every wf tail
    accepted_list = []          # count the accepted wf
    peak = findmax(wf.value)[2]
    base_start = wf.value[1]
    slope_t1 = tail_slope(wf, trunc(Int, (length(wf.value)-peak)*0.5))
    slope_t2 = tail_slope(wf, length(wf.value)-peak)
    slope_b  = base_slope(wf, 100)
    peak_index = risepoint(wf)
    baseline = wf.value[begin:peak_index]
    uplim, dolim = basestart(tier1_table.waveform)
    if (base_start < uplim && base_start > dolim) && (wf.value[1000] > 1000) && (peak < 2100 && peak > 1650) && (slope_t2 - slope_t1 < 0.18) && (base_start - mean(baseline) < 50) && (slope_b < 0.001 && slope_b > -0.001)
        append!(accepted_list, base_start)
        append!(std_list,std(baseline))
        slope_tau = tail_slope(wf, 1500)
        append!(tau_list, slope_tau)
    end
stats(sample, accepted_list, tau_list, std_list)
end

function selection_cut_2(wfs::ArrayOfRDWaveforms, sample)
    #main
    std_list = []               # store std for every wf
    tau_list = []               # store slope of every wf tail
    accepted_list = []          # count the accepted wf
    uplim, dolim = basestart(tier1_table.waveform)
    for i in 1:sample 
        wf = wfs[i]
        peak = findmax(wf.value)[2]
        slope_t1 = tail_slope(wf, trunc(Int, (length(wf.value)-peak)*0.5))
        slope_t2 = tail_slope(wf, length(wf.value)-peak)
        slope_b  = base_slope(wf, 800)
        base_start = wf.value[1]
        peak_index = risepoint(wf)
        baseline = wf.value[begin:peak_index]
        if (base_start < uplim && base_start > dolim) && (wf.value[1000] > 1000) && (peak < 2100 && peak > 1650) && (slope_t2 - slope_t1 < 0.18) && (base_start - mean(baseline) < 50) && (slope_b < 0.001 && slope_b > -0.001)
            append!(accepted_list, base_start)
            append!(std_list,std(baseline))
            slope_tau = tail_slope(wf, 1500)
            append!(tau_list, slope_tau)
            p = plot!(wf)
            display(p)
        end
    end
stats(sample, accepted_list, tau_list, std_list)
end
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
##
# wf = tier1_table.waveform[800]
# peak = findmax(wf.value)[2]
# plot(wf)
# plot!(peak_finder(wf, length(wf.value)-peak))
# plot!(peak_finder(wf, trunc(Int, (length(wf.value)-peak)*0.5)))
# plot!(extract_baseline(wf))
# ##
# function baseline_dist(wfs::ArrayOfRDWaveforms)       
#     #get the distribution of baseline starting point 
#     basestart_list = []
#     for i in wfs
#         append!(basestart_list, i.value[1])
#     end
#     return basestart_list
# end  
# #stephist(baseline_dist(tier1_table.waveform),  bins = 12000:1:13000, label = "Baseline start")
# ##
# #------------------------------------------------------------------------------
# #------------------------------------------------------------------------------
# ##
# #test
# #selection_cut(wf, 1)
# selection_cut_2(tier1_table.waveform, 10)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
##
##