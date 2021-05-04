"""
    mcpss_to_mcraw(mcpss, mctruth, sim_config_file)

Process simulated waveforms to account for DAQ and electronics effects,
resulting in a table that mimics raw data format.

mcpss: Table with simulated waveforms (output of mcstp_to_mcpss)
mctruth: Table with MC truth (currently used to add DAQ timestamp)
sim_config_file: simulation configuration json file
"""
function mcpss_to_mcraw(mcpss::Table, mctruth::Table, sim_config_file::AbstractString)
    # simulation config
    sim_config = LegendGeSim.load_config(sim_config_file)

    mcpss_to_mcraw(mcpss, mctruth, sim_config)
end   

"""
    mcpss_to_mcraw(mcpss, mctruth, sim_config)

Process simulated waveforms to account for DAQ and electronics effects,
resulting in a table that mimics raw data format.

mcpss: Table with simulated waveforms (output of mcstp_to_mcpss)
mctruth: Table with MC truth (currently used to add DAQ timestamp)
sim_config: PropDict object with simulation configuration
"""
function mcpss_to_mcraw(mcpss::Table, mctruth::Table, sim_config::PropDict)
    # simulation config
    daq = GenericDAQ(sim_config)
    preamp = PreAmp(sim_config)
    noise_model = NoiseModel(sim_config)

    mcpss_to_mcraw(mcpss, mctruth, daq, preamp, noise_model)
end    


"""
    mcpss_to_mcraw(mcpss, mctruth, daq, elec_chain, noise_model)

Process simulated waveforms to account for DAQ and electronics effects,
resulting in a table that mimics raw data format.

mcpss: Table with simulated waveforms (output of mcstp_to_mcpss)
mctruth: Table with MC truth (currently used to add DAQ timestamp)
daq: GenericDAQ object with DAQ parameters for the simulation
elec_chain: ElecChain object with parameters for the simulation of the electronics chain
noise_model: NoiseModel object for noise simulation
(deep simulaiton, sim based on data, baseline chunk slapping)

Output: Table

Work in progress
"""
function mcpss_to_mcraw(mcpss::Table, mctruth::Table, daq::GenericDAQ, elec_chain::ElecChain, noise_model::NoiseModel)

    result = process_waveforms(mcpss, daq, elec_chain, noise_model)
   
    # why am I doing this?
    wf_final = ArrayOfRDWaveforms(result.wf_array)
    wf_final = ArrayOfRDWaveforms((wf_final.time, VectorOfSimilarVectors(wf_final.value)))

    ##
    mcraw = Table(
        baseline = result.baseline,
        channel = mcpss.channel,
        energy = result.online_energy,
        ievt = mcpss.ievt,
        # why am i using baseline here? could be any column?
        numtraces = ones(length(result.baseline)), # number of triggered detectors (1 for HADES)
        packet_id = zeros(length(result.baseline)), # means to packet losses
        timestamp = getindex.(mctruth.thit, 1), # frist MC truth hit time of each event?
        tracelist = VectorOfVectors([[1] for idx in 1:length(result.baseline)]), # lists of ADCs that triggered, 1 for HADES all the time
        waveform = wf_final,
        wf_max = maximum.(wf_final.value),
        # I was told that std is std of the whole waveform, not just the baseline
        wf_std = std.(wf_final.value)
    )

    # probably not needed because one can look this up in the simulation config file
    # electruth = Table(
    #     # fill in noise and DAQ parameters
    # )

    mcraw

end


function process_waveforms(mcpss::Table, daq::GenericDAQ, elec_chain::ElecChain, noise_model::NoiseSim)
    # Create arrays to be filled with results, and online energy
    n_waveforms = size(mcpss.waveform,1)
    result = Table(
        wf_array = Array{RDWaveform}(undef, n_waveforms),
        online_energy = Array{typeof(1.0*energy_unit)}(undef, n_waveforms),
        baseline = Array{T}(undef, n_waveforms),
        baseline_rms = Array{T}(undef, n_waveforms)
    )

    @info "Processing waveforms..."
    for i in 1:n_waveforms
        # @showprogress 1 "Processing..." for i in 1:n_waveforms
        if(i % 500 == 0) println("$i / $n_waveforms") end
        # println("$i / $n_waveforms")

        # plot_wf(mcpss.waveform[i],1)

        wf = add_tail_and_baseline(mcpss.waveform[i], daq)
        # plot_wf(wf_array[i],2)

        
        wf = differentiate(wf)
        # plot_wf(wf_array[i],3)


        wf = simulate_elec(wf, elec_chain)
        # plot_wf(wf_array[i],4)

        # works for preamp noise, but not data noise, since that's in DAQ units...
        wf = simulate_noise(wf, noise_model)
        # plot_wf(wf_array[i],5)    

        wf = simulate_daq(wf, daq)
        # plot_wf(wf_array[i],6)

        # if online energy is zero, means didn't trigger
        # trigger_threshold = noise_σ * 3 * daq.gain # noise implemented before DAQ gain
        result.wf_array[i], result.online_energy[i] = trigger(wf, daq, noise_model.noise_σ*daq.gain*3)
        # plot_wf(wf_array[i],7)


       result.baseline[i], result.baseline_rms[i] = mean_and_std(result.wf_array[i].value[1:daq.baseline_length])

    end
    result
end


function process_waveforms(mcpss::Table, daq::GenericDAQ, elec_chain::ElecChain, noise_model::NoiseData)
    # Create arrays to be filled with results, and online energy
    n_waveforms = size(mcpss.waveform,1)
    result = Table(
        wf_array = Array{RDWaveform}(undef, n_waveforms),
        online_energy = Array{typeof(1.0*energy_unit)}(undef, n_waveforms),
        baseline = Array{T}(undef, n_waveforms),
        baseline_rms = Array{T}(undef, n_waveforms)
    )

    # pick a random baseline from catalog
    baselines_table = baseline_catalog(noise_model.baseline_catalog)
    baseline = rand(Tables.getcolumn(Tables.columns(baselines_table), :waveform))
    Δt = step(baseline.time)
    if Δt > daq.Δt
        println("Time step in baseline sample ($Δt) is larger than DAQ Δt ($(daq.Δt)) based on given sampling rate!")
        return
    end

    @info "Processing waveforms..."
    for i in 1:n_waveforms
        if(i % 500 == 0) println("$i / $n_waveforms") end

        wf = add_tail_and_baseline(mcpss.waveform[i], daq)
        
        wf = differentiate(wf)

        wf = simulate_elec(wf, elec_chain)
        # plot_wf(wf_array[i],4)

        wf = simulate_daq(wf, daq, baseline)
        # plot_wf(wf_array[i],6)

        # if online energy is zero, means didn't trigger
        result.wf_array[i], result.online_energy[i] = trigger(wf, daq, std(baseline.value)*3)
        # plot_wf(wf_array[i],7)

       result.baseline[i], result.baseline_rms[i] = mean_and_std(result.wf_array[i].value[1:daq.baseline_length])
    end

    result
end


"""
    differentiate(wf)

Differentiate a waveform using Biquad filter
(see function dspjl_differentiator_filter)

wf: RDWaveform
Output: RDWaveform
"""
function differentiate(wf::RDWaveform)
    diff_biquad_filter = dspjl_differentiator_filter(T(1)) # We want a gain of 1 here
    filter_output = filt(diff_biquad_filter, wf.value)
    # !!! We have to implement a method do parse a RDWaveform to `filt`
    RDWaveform(wf.time, filter_output)
end



##

"""
    add_tail_and_baseline(wf, daq)

Extend tail and add baseline based on DAQ number of samples and baseline length
To be changed to take ArrayOfRDWaveforms rather than a single RDWaveform

wf: RDWaveform object 
daq: DAW object
"""
function add_tail_and_baseline(wf::RDWaveform, daq::DAQ)
# function add_tail_and_baseline(wfs::RadiationDetectorSignals.ArrayOfRDWaveforms, daq::DAQ)
    factor::Integer = uconvert(u"ns", daq.Δt) / uconvert(u"ns", step(wf.time))
    SolidStateDetectors.add_baseline_and_extend_tail(wf, daq.baseline_length*factor*2, daq.nsamples*factor*2)
    # ArrayOfRDWaveforms(SolidStateDetectors.add_baseline_and_extend_tail.(wfs, daq.baseline_length*factor, daq.nsamples*factor))
    # SolidStateDetectors.add_baseline_and_extend_tail.(wfs, daq.baseline_length*factor, daq.nsamples*factor)

end




"""
    trigger(wf, daq, noise_model)

Simulate DAQ trigger.
Returns the waveform stored by the DAQ and the resulting online energy

wf: RDWaveform
daq: GenericDAQ object
noise_model: NoiseModel object

Output: RDWaveform, float
"""
# function trigger(wf::RDWaveform, daq::GenericDAQ, noise_σ::Real)
function trigger(wf::RDWaveform, daq::GenericDAQ, trigger_threshold::Real)
    trigger_window_lengths = (250,250,250)
    trigger_window_length = sum(trigger_window_lengths)
    # trigger_threshold = noise_σ * 3 * daq.gain # noise implemented before DAQ gain

    online_filter_output = zeros(T, length(wf.value) - trigger_window_length)
    t0_idx::Int = 0
    trig = false

    for i in eachindex(online_filter_output)
        online_filter_output[i], trig = daq_online_filter(wf.value, i, trigger_window_lengths, trigger_threshold)
        if trig && t0_idx == 0
            t0_idx = i+trigger_window_lengths[1]+trigger_window_lengths[2]
        end
    end

    ts = range(T(0)u"ns", step = daq.Δt, length = daq.nsamples)
    # in case it didn't trigger
    if(t0_idx == 0)
        stored_waveform = RDWaveform(ts, wf.value[1:daq.nsamples]) # just to return something
        online_energy = 0u"keV" # flag meaning didn't trigger
    else

        online_energy = uconvert(u"keV", maximum(online_filter_output) * germanium_ionization_energy / daq.gain)

        # plot_wf = plot(wf, xlims=(0,20000))
        # vline!([ustrip(wf.time[t0_idx]),
        #         ustrip(wf.time[t0_idx - trigger_window_lengths[2]]),
        #         ustrip(wf.time[t0_idx - trigger_window_lengths[2] - trigger_window_lengths[1]]),
        #         ustrip(wf.time[t0_idx+trigger_window_lengths[3]])])
        # png(plot_wf, "wf_7a.png")

        iStart = t0_idx - daq.baseline_length

        # wf1 = RDWaveform(wf.time[t_:end], wf.value[iStart:end])
        # plot!(wf1, linestype=:dot, linecolor=:red)
        # iStart = t0_idx-daq.baseline_length*Int(daq.Δt/step(wv.time))

        stored_waveform = RDWaveform(ts, wf.value[iStart:iStart+daq.nsamples-1])

    end

    stored_waveform, online_energy
end

# function plot_wf(wf, idx)
#     plt_wf = plot(wf, xlims=(0,20000))
#     png(plt_wf, "wf_$idx.png")    
# end