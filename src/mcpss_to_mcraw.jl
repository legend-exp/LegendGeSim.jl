keV_unit = 1.0u"keV"; keV = typeof(keV_unit);
MHz_unit = 1.0u"MHz"; MHz = typeof(MHz_unit);
ns_unit = 1.0u"ns"; ns = typeof(ns_unit);
μs_unit = 1.0u"μs"; μs = typeof(μs_unit);


## Fix basic parameters

T = Float32

germanium_ionization_energy = T(2.95)u"eV"

"""
A shelf to put all the DAQ parameters in
later to be filled from a configuration file
"""
@with_kw struct DAQ
    "How many samples DAQ stores (final wf length). Has to be smaller than total wf length (defined in mcraw_to_mcpss.jl)"
    nsamples::Int = 4000; # samples per waveform

    "Length of DAQ baseline"
    baseline_length::Int = 400;

    "DAQ sampling rate. Needed to calculate Δt"
    sampling_rate::MHz = 250u"MHz"

    "Inverse of sampling rate. Needed for DAQ trigger"
    Δt::ns = uconvert(u"ns", inv(sampling_rate));

    "DAQ type. Needed to make sure the wf values are..."
    daq_type = UInt16

    "Maximum energy DAQ can register?"
    max_e::keV = 10_000u"keV" # == typemax(UInt16)

    "DAQ offset"
    offset::keV = 2_000u"keV";

    "What is this parameter?"
    c = typemax(daq_type) / ((max_e+offset)/uconvert(u"keV", germanium_ionization_energy))
end

daq = DAQ()


"""
A shelf to put all the preamplifier parameters in
later to be filled from a configuration file
"""
@with_kw mutable struct PreAmp
    "PreAmp exp decay time"
    τ_decay::μs = (50)*u"μs"
    
    "PreAmp rise time"
    τ_rise::ns = T(15)*u"ns"
end

pa = PreAmp()

# sigma for electronic noise
const noise_σ = uconvert(u"eV", T(3)u"keV") / germanium_ionization_energy

##

"""
    mcpss_to_mcraw(mcpss, mctruth)

Process simulated waveforms to account for DAQ and electronics effects,
resulting in a table that mimics raw data format.

mcpss: Table with simulated waveforms (output of mcstp_to_mcpss)
mctruth: Table with MC truth (currently used to add DAQ timestamp)

Output: Table
"""
function mcpss_to_mcraw(mcpss::Table, mctruth::Table)
    ### Create arrays to be filled with results, and online energy
    idx_end = size(mcpss.waveform,1)
    wf_array = Array{RDWaveform}(undef, idx_end)
    online_energy = Array{keV}(undef, idx_end)
    baseline = Array{T}(undef, idx_end)
    baseline_rms = Array{T}(undef, idx_end)

    @info "Processing waveforms..."
    ### loop over each wf and process it
   for i in 1:idx_end
    # @showprogress 1 "Processing..." for i in 1:idx_end
       if(i % 500 == 0) println("$i / $idx_end") end
#        println("$i / $idx_end")

        ### Differentiate
        wf = differentiate_wf(mcpss.waveform[i])
        wf_array[i] = wf

        ### 1. PreAmp Simulation

        #### 1.1. Simulate a charge sensitive amplifier (`CSA`)
        wf_array[i] = simulate_csa(wf_array[i])

        #### 1.2. Noise
        wf_array[i] = simulate_noise(wf_array[i])

        ### 2. DAQ Simulation

        # Now we want to simulate, what the DAQ records
        # The waveform will go rough the DAQ buffer on which the DAQ filter runs
        # in order to trigger in case of an event and calculate some parameters
        # like online energy
        # In order to do so. The waveform should be longer the than the DAQ Buffer (here, 5000 samples)
        # We took care of this at the beginning

        #### 2.1. DAQ units and baseline
        wf_array[i] = daq_baseline(wf_array[i])

        #### 2.2. Trigger method
        # if online energy is zero, means didn't trigger
        wf_array[i], online_energy[i] = daq_trigger(wf_array[i])

        baseline[i], baseline_rms[i] = mean_and_std(wf_array[i].value[1:daq.baseline_length])

    end

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
    simulate_csa(wf)

Simulate a charge sensitive amplifier (`CSA`)
Here, the parameters τ_rise and τ_decay have to be given in units of samples,
because the `filt`-function does not know the Δt between the samples.

wf: RDWaveform
Output: RDWaveform
"""
function simulate_csa(wf::SolidStateDetectors.RDWaveform)

    csa_filter = dspjl_simple_csa_response_filter(
        pa.τ_rise / step(wf.time),
        uconvert(u"ns", pa.τ_decay) / step(wf.time))

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


##
"""
    simulate_noise(wf)

Simulate electronic noise
Currently simple Gaussian noise. Parameters defined on top of the script
(currently dummy values).
TODO: parameters read from config file

wf: RDWaveform
Output: RDWaveform
"""
function simulate_noise(wf::SolidStateDetectors.RDWaveform)
    # I am no expert here. I don't know at which point one should introduce noise.
    # Also, different noise could be added at different stages. This really depends on the electronics.
    # I will just add some Gaussian Noise (σ of 3 keV defined on top)
    # lets generate 1000 random samples from this normal distribution
    gaussian_noise_dist = Normal(T(0), T(noise_σ)) #  Normal() -> Distributions.jjl
    samples = rand(gaussian_noise_dist, 1000)

    # Now, lets add this Gaussian noise to other waveform (here, after the filters (but might be also added before))
    wf_noise = RDWaveform(wf.time, wf.value .+ rand!(gaussian_noise_dist, similar(wf.value)))
    wf_noise
end

##
"""
    daq_baseline(wf)

Add DAQ baseline

wf: RDWaveform
Output:
"""
function daq_baseline(wf::SolidStateDetectors.RDWaveform)
    o = daq.c * uconvert(u"eV", daq.offset) / germanium_ionization_energy

    # invert the pulse if needed
    sign = wf.value[end] < 0 ? -1 : 1

    # daq_buffer_wv = RDWaveform(wf.time, daq.daq_type.(round.(sign * wf.value .* daq.c .+ o, digits = 0)))
    daq_buffer_wv = RDWaveform(wf.time, daq.daq_type.(round.(sign / ustrip(germanium_ionization_energy) * wf.value .* daq.c .+ o, digits = 0)))
    daq_buffer_wv
end


"""
    daq_trigger(wf)

Simulate DAQ trigger. Returns a waveform with trigger and resulting online energy

wf: RDWaveform
Output: RDwaveform, float
"""
function daq_trigger(wf::SolidStateDetectors.RDWaveform)
    daq_trigger_window_lengths = (250,250,250)
    daq_trigger_window_length = sum(daq_trigger_window_lengths)
    daq_trigger_threshold = noise_σ * 10 * daq.c

    online_filter_output = zeros(T, length(wf.value) - daq_trigger_window_length)
    t0_idx = 0
    trig = false

    for i in eachindex(online_filter_output)
        online_filter_output[i], trig = daq_online_filter(wf.value, i-1, daq_trigger_window_lengths, daq_trigger_threshold)
        if trig && t0_idx == 0
            t0_idx = i
        end
    end

    ts = range(T(0)u"ns", step = daq.Δt, length = daq.nsamples)
    # in case it didn't trigger
    if(t0_idx == 0)
        stored_waveform = RDWaveform(ts, wf.value[1:daq.nsamples]) # just to return something
        online_energy = 0u"keV" # flag meaning didn't trigger
    else
        online_energy = uconvert(u"keV", maximum(online_filter_output) * germanium_ionization_energy / daq.c)
        iStart = t0_idx-daq.baseline_length
        stored_waveform = RDWaveform(ts, wf.value[iStart:iStart+daq.nsamples-1]);
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
