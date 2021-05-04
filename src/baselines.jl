function extend_baseline(baseline::RDWaveform, wf::RDWaveform)
    # resample baseline
    baseline_sampled = baseline.value[begin:Int(step(wf.time)/step(baseline.time)):end]

    # create extended/shrinked baseline
    gaussian_noise_dist = Normal(T(0), T(std(baseline_sampled))) #  Normal() -> Distributions.jjl
    values = ones(length(wf.value)).*mean(baseline_sampled)
    values = values .+ rand!(gaussian_noise_dist, similar(values))
    # Now, lets add this Gaussian noise to other waveform (here, after the filters (but might be also added before))
    RDWaveform(wf.time, values)
end


function baseline_catalog(raw_filename::AbstractString)
    base_filename = joinpath("cache", "baselines_"*basename(raw_filename))

    if isfile(base_filename)
        @info "Selecting baseline samples from $base_filename"
        baseline_table = HDF5.h5open(base_filename) do input Table(waveform = LegendDataTypes.readdata(input, "raw/waveform")) end
    else
        @info "Extracting baseline samples from $raw_filename"    
        raw_table = read_raw(raw_filename, "raw")
        baseline_table = baseline_catalog(raw_table)
        HDF5.h5open(base_filename, "w") do f LegendDataTypes.writedata(f, "raw", baseline_table) end
        @info "Baselines saved to $base_filename"
    end

    baseline_table
end


function baseline_catalog(raw_table::Table)
    waveforms = raw_table.waveform
    # Tns = typeof(1.0*ns_unit)
    # times = Vector{Vector{Tns}}()
    # values = Vector{Vector{T}}()
    baselines = Vector{RDWaveform}()

    base_uplim, base_lolim = basestart(waveforms)

    for wf in waveforms
        if(selection_cut(wf, base_uplim, base_lolim))
            push!(baselines, extract_baseline(wf))
            # time, value = extract_baseline(wf)
            # push!(times, time)
            # push!(values, value)
        end
    end

    # ArrayOfRDWaveforms((VectorOfVectors(times), VectorOfVectors(waveforms)))
    Table(waveform = ArrayOfRDWaveforms(baselines))
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


function selection_cut(wf::RDWaveform, base_uplim::Real, base_lolim::Real)
    # baseline start cut
    base_start = wf.value[1]
    peak_index = risepoint(wf)
    baseline = wf.value[begin:peak_index]
    cut_base::Bool = (base_start < base_uplim) && (base_start > base_lolim) && (base_start - mean(baseline) < 50)

    # wveform value cut 
    cut_value = wf.value[1000] > 1000

    # peak cut 
    peak = findmax(wf.value)[2]
    cut_peak = peak < 2100 && peak > 1650

    # slope cut 
    slope_t1 = tail_slope(wf, trunc(Int, (length(wf.value)-peak)*0.5))
    slope_t2 = tail_slope(wf, length(wf.value)-peak)
    cut_slope = slope_t2 - slope_t1 < 0.18

    # baseline slope cut
    slope_b  = base_slope(wf, 100)
    cut_bslope = slope_b < 0.001 && slope_b > -0.001

    cut_base && cut_value && cut_peak && cut_slope && cut_bslope
end


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


function tail_slope(wf::RDWaveform, n)
    #slope of the tail, n = 500, n= 1500, n=2000
    peak = findmax(wf.value)[2]
    charge = wf.value[peak:peak+n]
    time = wf.time[peak:peak+n]
    _, slope = linear_fit(ustrip(time), ustrip(charge))
    slope
end



function base_slope(wf::RDWaveform, n)
    #slope of the tail  n = 1800
    peak = findmax(wf.value)[2]
    charge = wf.value[begin:peak-n]
    time = wf.time[begin:peak-n]
    _, slope = linear_fit(ustrip(time), ustrip(charge))
    slope
end


function extract_baseline(wf::RDWaveform)
    #extract baseline
    baseline_index = risepoint(wf)
    # wf.time[begin:baseline_index], wf.value[begin:baseline_index]
    RDWaveform(wf.time[begin:baseline_index], wf.value[begin:baseline_index])
end

# # after we read this function, we obtain a table
# # one of its columns is "waveform" and it contains RDWaveform objects
# filename = joinpath(tier1_path, tier1_file)
# tier1_table = read_raw(filename, "raw")
# ##
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------




# function peak_finder(wf::RDWaveform,n)
#     peak = findmax(wf.value)[2]
#     return RDWaveform(wf.time[peak:peak+n], wf.value[peak:peak+n])
# end








# function stats(sample, accepted_list, tau_list, stdlist)
#     #stats
# println("No. of waveforms = ", sample) 
# println("Accepted = ", length(accepted_list), "  ----------  ", (length(accepted_list)/sample)*100 , " %") 
# println("mean std = ", mean(stdlist))                                  #means std
# println("max tau = ", -1/maximum(tau_list)*16/1000 , " microsecond")    #max tau
# println("min tau = ", -1/minimum(tau_list)*16/1000 , " microsecond")    #min tau
# end



# function selection_cut_2(wfs::ArrayOfRDWaveforms, sample)
#     #main
#     std_list = []               # store std for every wf
#     tau_list = []               # store slope of every wf tail
#     accepted_list = []          # count the accepted wf
#     uplim, dolim = basestart(tier1_table.waveform)
#     for i in 1:sample 
#         wf = wfs[i]
#         peak = findmax(wf.value)[2]
#         slope_t1 = tail_slope(wf, trunc(Int, (length(wf.value)-peak)*0.5))
#         slope_t2 = tail_slope(wf, length(wf.value)-peak)
#         slope_b  = base_slope(wf, 800)
#         base_start = wf.value[1]
#         peak_index = risepoint(wf)
#         baseline = wf.value[begin:peak_index]
#         if (base_start < uplim && base_start > dolim) && (wf.value[1000] > 1000) && (peak < 2100 && peak > 1650) && (slope_t2 - slope_t1 < 0.18) && (base_start - mean(baseline) < 50) && (slope_b < 0.001 && slope_b > -0.001)
#             append!(accepted_list, base_start)
#             append!(std_list,std(baseline))
#             slope_tau = tail_slope(wf, 1500)
#             append!(tau_list, slope_tau)
#             p = plot!(wf)
#             display(p)
#         end
#     end
# stats(sample, accepted_list, tau_list, std_list)
# end
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