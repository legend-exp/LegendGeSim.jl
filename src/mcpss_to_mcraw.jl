"""
    mcpss_to_mcraw(mcpss, mctruth)

Process simulated waveforms to account for DAQ and electronics effects,
resulting in a table that mimics raw data format.

mcpss: Table with simulated waveforms (output of mcstp_to_mcpss)
mctruth: Table with MC truth (currently used to add DAQ timestamp)

Output: Table
"""
function mcpss_to_mcraw(mcpss::Table, mctruth::Table, sim_config::PropDict)
    # function mcpss_to_mcraw(mcpss::Table, mctruth::Table, sim_conf_file::AbstractString)
    ### Create arrays to be filled with results, and online energy
    n_waveforms = size(mcpss.waveform,1)
    wf_array = Array{RDWaveform}(undef, n_waveforms)
    online_energy = Array{typeof(1.0*energy_unit)}(undef, n_waveforms)
    baseline = Array{T}(undef, n_waveforms)
    baseline_rms = Array{T}(undef, n_waveforms)

    ## electronics configuration
    # @info "Elec config taken from: $sim_conf_file"
    # sim_conf = PropDicts.read(PropDict, sim_conf_file)
    daq = construct_GenericDAQ(sim_config)
    preamp = construct_PreAmp(sim_config)

    @info "Processing waveforms..."
    ### loop over each wf and process it
   for i in 1:n_waveforms
    # @showprogress 1 "Processing..." for i in 1:n_waveforms
       if(i % 500 == 0) println("$i / $n_waveforms") end
#        println("$i / $n_waveforms")

        ### Differentiate
        wf_array[i] = differentiate_wf(mcpss.waveform[i])

        ### Extend tail and add baseline
        wf_array[i] = tail_and_baseline(wf_array[i], daq)

        ### 1. PreAmp Simulation

        #### 1.1. Simulate a charge sensitive amplifier (`CSA`)
        wf_array[i] = simulate_csa(wf_array[i], preamp)

        #### 1.2. PreAmp/Data noise
        noise_sigma = sim_config.noise_data == 0 ? preamp.noise_σ : sim_config.noise_data
        # println("Noise: $noise_sigma")
        wf_array[i] = simulate_noise(wf_array[i], noise_sigma)


        ### 2. DAQ Simulation

        # Now we want to simulate, what the DAQ records
        # The waveform will go rough the DAQ buffer on which the DAQ filter runs
        # in order to trigger in case of an event and calculate some parameters
        # like online energy
        # In order to do so. The waveform should be longer the than the DAQ Buffer (here, 5000 samples)
        # We took care of this at the beginning

        #### 2.1 Gain & offset

        wf_array[i] = amplify_and_offset(wf_array[i], daq)

        #### 2.2. Resample & digitize

        wf_array[i] = resample_and_digitize(wf_array[i], daq)

        #### 3.2. Trigger method
        # if online energy is zero, means didn't trigger
        wf_array[i], online_energy[i] = trigger(wf_array[i], daq)

        baseline[i], baseline_rms[i] = mean_and_std(wf_array[i].value[1:daq.baseline_length])

    end

    # why am I doing this?
    wf_final = ArrayOfRDWaveforms(wf_array)
    wf_final = ArrayOfRDWaveforms((wf_final.time, VectorOfSimilarVectors(wf_final.value)))


    ##
    mcraw = Table(
        baseline = baseline,
        channel = mcpss.channel,
        energy = online_energy,
        ievt = mcpss.ievt,
        numtraces = ones(length(baseline)), # number of triggered detectors (1 for HADES)
        packet_id = zeros(length(baseline)), # means to packet losses
        timestamp = getindex.(mctruth.thit, 1), # frist MC truth hit time of each event?
        tracelist = VectorOfVectors([[1] for idx in 1:length(baseline)]), # lists of ADCs that triggered, 1 for HADES all the time
        waveform = wf_final,
        wf_max = maximum.(wf_final.value), # ?
        # wf_std = Vector([ std(wf_final[idx].value[1:daq.baseline_length] for idx in 1:length(wf_final))])
        wf_std = std.(wf_final.value) # ?
        # wf_std = std.(wf_final.value[1:daq.baseline_length]) # ?
    )

    # electruth = Table(
    #     # fill in noise and DAQ parameters
    # )

    mcraw


end


"""
    differentiate_wf(wf)

Differentiate a waveform using Biquad filter
(see function dspjl_differentiator_filter)

wf: RDWaveform
Output: RDWaveform
"""
function differentiate_wf(wf::SolidStateDetectors.RDWaveform)
    diff_biquad_filter = dspjl_differentiator_filter(T(1)) # We want a gain of 1 here
    filter_output = filt(diff_biquad_filter, wf.value)
    # !!! We have to implement a method do parse a RDWaveform to `filt`
    wf_diff = RDWaveform(wf.time, filter_output)
    wf_diff
end

"""
    dspjl_differentiator_filter(gain::Real)

Differentiator filter (later to be part of LegendDSP)

The waveform is in units of induced charge (integrated form).
But we want to transform it into its differential form (current).
We will use a BiQuad filter for this: https://en.wikipedia.org/wiki/Digital_biquad_filter

gain: ?
Output: ?
"""
function dspjl_differentiator_filter(gain::Real)

    T = float(typeof(gain))
    Biquad(T(gain), T(-gain), T(0), T(0), T(0))
end


##

"""
    tail_and_baseline(wf, daq)
    Extend tail and add baseline based on DAQ number of samples and baseline length
"""
function tail_and_baseline(wf::SolidStateDetectors.RDWaveform, daq::DAQ)
    # WRONG, should not use directly daq params, but convert accounting for 1ns/4ns
    wf_tail = SolidStateDetectors.add_baseline_and_extend_tail(wf, daq.baseline_length, daq.nsamples)
    wf_tail
end

##

"""
    simulate_csa(wf)

Simulate a charge sensitive amplifier (`CSA`)
Here, the parameters τ_rise and τ_decay have to be given in units of samples,
because the `filt`-function does not know the Δt between the samples.

wf: RDWaveform
Output: RDWaveform
"""
function simulate_csa(wf::SolidStateDetectors.RDWaveform, preamp::PreAmp)

    csa_filter = dspjl_simple_csa_response_filter(
        preamp.τ_rise / step(wf.time),
        uconvert(u"ns", preamp.τ_decay) / step(wf.time))

    pa_wf = RDWaveform(wf.time, filt(csa_filter, wf.value))
    pa_wf
end


"""
    dspjl_simple_csa_response_filter(τ_rise, τ_decay, gain)

Simulate CSA response using the RC and integrator filters
τ_rise: ?
τ_decay: ?
gain: ?

Output: ?
"""
function dspjl_simple_csa_response_filter(τ_rise::Real, τ_decay::Real, gain::Real = one(τ_rise))
    # TODO: Use a single biquad filter
    T = float(promote_type(promote_type(typeof(τ_rise), typeof(τ_decay)), typeof(gain)))
    dspjl_rc_filter(T(τ_rise)) * dspjl_integrator_cr_filter(T(gain), T(τ_decay))
end


"""
    dspjl_rc_filter(RC)

An `RC` filter made with BiQuad filter

Differentiate a waveform using Biquad filter
(see function dspjl_differentiator_filter)

wf: RDWaveform
Output: RDWaveform
"""
function dspjl_rc_filter(RC::Real)
    T = float(typeof(RC))
    α = 1 / (1 + RC)
    Biquad(T(α), T(0), T(0), T(α - 1), T(0))
end

"""
    dspjl_integrator_cr_filter(gain, RC)

An `integrator` filter (the inverser of the `dspjl_differentiator_filter`)

gain: ?
RC: ?

Output: ?
"""
function dspjl_integrator_cr_filter(gain::Real, RC::Real)
    T = float(promote_type(typeof(gain), typeof(RC)))
    α = 1 / (1 + RC)
    Biquad(T(gain), T(-α), T(0), T(α - 1), T(0))
end


"""
    simulate_noise(wf, noise_σ)

Simulate noise based on sigma extracted from data.
Should be applied after conversion to DAQ units.
Currently simple Gaussian noise.

wf: RDWaveform
noise_σ: value in DAQ units
Output: RDWaveform
"""
function simulate_noise(wf::SolidStateDetectors.RDWaveform, noise_σ::Real)
    # I am no expert here. I don't know at which point one should introduce noise.
    # Also, different noise could be added at different stages. This really depends on the electronics.
    # I will just add some Gaussian Noise (σ of 3 keV defined on top)
    # lets generate 1000 random samples from this normal distribution
    # println("Noise: $noise_σ")
    gaussian_noise_dist = Normal(T(0), T(noise_σ)) #  Normal() -> Distributions.jjl
    samples = rand(gaussian_noise_dist, 1000)

    # Now, lets add this Gaussian noise to other waveform (here, after the filters (but might be also added before))
    wf_noise = RDWaveform(wf.time, wf.value .+ rand!(gaussian_noise_dist, similar(wf.value)))
    wf_noise
end


##
"""
    amplify_and_offset(wf, daq)

Gain & offset (currently in keV)

wf: RDWaveform
daq: DAQ object
Output: RDWaveform
"""
function amplify_and_offset(wf::SolidStateDetectors.RDWaveform, daq::DAQ)
    o = daq.gain * uconvert(u"eV", daq.offset) / germanium_ionization_energy

    # invert the pulse if needed
    sign = wf.value[end] < 0 ? -1 : 1

    # daq_buffer_wv = RDWaveform(wf.time, daq.daq_type.(round.(sign * wf.value .* daq.gain .+ o, digits = 0)))
    # daq_buffer_wv = RDWaveform(wf.time, UInt16.(round.(sign / ustrip(germanium_ionization_energy) * wf.value .* daq.gain .+ o, digits = 0)))
    daq_buffer_wv = RDWaveform(wf.time, sign / ustrip(germanium_ionization_energy) * wf.value .* daq.gain .+ o)
    daq_buffer_wv
end

##

"""
    resample_and_digitize(wf, daq)

SSD samples are 1ns apart -> resample based on DAQ Δt
Afterwards digitize (currently generic DAQ has 16 (14?) bits)

wf: RDWaveform
daq: GenericDAQ object
Output: RDWaveform
"""
function resample_and_digitize(wf::SolidStateDetectors.RDWaveform, daq::GenericDAQ)
    ts = range(T(0)u"ns", step = daq.Δt, length = daq.nsamples)
    # resample
    wf_samp = RDWaveform(ts, wf.value[range(0, step = Int(daq.Δt/step(wf.time)), length = daq.nsamples)])
    # digitize
    wf_dig = RDWaveform(wf_dig.time, UInt16.(round.(wf_samp.value, digits = 0)))
    wf_dig

end


"""
    trigger(wf)

Simulate DAQ trigger. Returns a waveform with trigger and resulting online energy

wf: RDWaveform
Output: RDwaveform, float
"""
function trigger(wf::SolidStateDetectors.RDWaveform, daq::DAQ)
    trigger_window_lengths = (250,250,250)
    trigger_window_length = sum(daq_trigger_window_lengths)

    online_filter_output = zeros(T, length(wf.value) - trigger_window_length)
    t0_idx = 0
    trig = false

    for i in eachindex(online_filter_output)
        online_filter_output[i], trig = daq_online_filter(wf.value, i-1, trigger_window_lengths, daq.trigger_threshold)
        if trig && t0_idx == 0
            t0_idx = i
        end
    end

    # ts = range(T(0)u"ns", step = daq.Δt, length = daq.nsamples)
    # in case it didn't trigger
    if(t0_idx == 0)
        stored_waveform = RDWaveform(ts, wf.value[1:daq.nsamples]) # just to return something
        online_energy = 0u"keV" # flag meaning didn't trigger
    else
        online_energy = uconvert(u"keV", maximum(online_filter_output) * germanium_ionization_energy / daq.gain)
        # iStart = t0_idx-daq.baseline_length
        iStart = t0_idx-daq.baseline_length*Int(daq.Δt/step(wv.time))

        # stored_waveform = RDWaveform(ts, wf.value[iStart:iStart+daq.nsamples-1]);
        stored_waveform = RDWaveform(ts, wf.value[range(iStart, step = Int(daq.Δt/step(wf.time)), length = daq.nsamples)]);

    end

    stored_waveform, online_energy
end

"""
    daq_online_filter(values, offset, window_lengths, threshold)

Used to simulate DAQ output (?)
(see function daq_trigger)

values: waveform values
offset: ?
window_lengths: ?
threshold: ?

Output: ?
"""
function daq_online_filter(values::AbstractVector, offset::Int, window_lengths::NTuple{3, Int}, threshold)
    wl = window_lengths
    r1 = offset+wl[1]+wl[2]:offset+wl[1]+wl[2]+wl[3]
    r2 = offset+wl[1]+wl[2]+wl[3]:offset+sum(window_lengths)
    r = mean(values[r2]) - mean(values[r1])
    r, r >= threshold
end




#### Should not be here:
##

"""
    read_mcpss(filename)

Helper function to read simulation output from an HDF5 file with mcpss format
(to be rewritten as a mcstp_to_mcpss() function that reads from storage)

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
(to be rewritten as a mcstp_to_mcpss() function that reads from storage)

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
##