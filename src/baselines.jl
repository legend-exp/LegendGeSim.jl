"""
    extend_baseline(baseline, wf)

RDWaveform, RDWaveform -> RDWaveform

Take a given <baseline> and extend it to match the length (and sampling)        
    of the given waveform <wf>
"""
function extend_baseline(baseline::RDWaveform, wf::RDWaveform)
    # resample baseline to be the same as the waveform
    baseline_sampled = baseline.value[begin:Int(step(wf.time)/step(baseline.time)):end]

    ## create extended/shrinked baseline
    # offset
    values = ones(length(wf.value)).*mean(baseline_sampled)
    # noise
    gaussian_noise_dist = Normal(T(0), T(std(baseline_sampled)))
    values = values .+ rand!(gaussian_noise_dist, similar(values))

    RDWaveform(wf.time, values)
end


"""
    baseline_catalog(raw_filename)

AbstractString -> Table 

Look up stored table of baselines corresponding to given raw data <raw_filename>.
If does not exist, construct such table.
"""
function baseline_catalog(raw_filename::AbstractString)
    base_filename = joinpath("cache", "baselines_"*basename(raw_filename))

    if isfile(base_filename)
        @info "Selecting baseline samples from $base_filename"
        baseline_table = h5open(base_filename) do input Table(waveform = LegendDataTypes.readdata(input, "raw/waveform")) end
    else
        @info "Extracting baseline samples from $raw_filename"    
        raw_table = read_raw(raw_filename, "raw")
        baseline_table = baseline_catalog(raw_table)
        # cache for later
        h5open(base_filename, "w") do f writedata(f, "raw", baseline_table) end
        @info "Baselines saved to $base_filename"
    end

    baseline_table
end

"""
    baseline_catalog(raw_table)

Table -> Table 

Construct table of baselines extracted from the waveforms
    contained in the given raw data table <raw_table>
"""
function baseline_catalog(raw_table::Table)
    waveforms = raw_table.waveform
    baselines = Vector{RDWaveform}()

    # upper and lower limit on initial offset value
    offset_uplim, offset_lolim = basestart(waveforms)

    for wf in waveforms
        # if the waveform passes the selection cuts for pileup removal...
        if(selection_cut(wf, offset_uplim, offset_lolim))
            # ...extract its baseline and collect it in the list
            push!(baselines, extract_baseline(wf))
        end
    end

    Table(waveform = ArrayOfRDWaveforms(baselines))
end    


function basestart(wfs::ArrayOfRDWaveforms)       
    # get the distribution of offset values at start 
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
    RDWaveform(wf.time[begin:baseline_index], wf.value[begin:baseline_index])
end