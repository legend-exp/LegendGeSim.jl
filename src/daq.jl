"""
Abstract DAQ type for hierarchy and multiple dispatch
"""
abstract type DAQ
end


"""
A shelf to put all the DAQ parameters in
"""
@with_kw mutable struct GenericDAQ <: DAQ
    "How many samples DAQ stores (final wf length). Has to be smaller than total wf length (defined in mcraw_to_mcpss.jl)"
    nsamples::Int = 4000; # samples per waveform

    "Length of DAQ baseline in number of samples"
    baseline_length::Int = 400;

    "Inverse of sampling rate. Needed for DAQ trigger"
    Δt::typeof(1.0*ns_unit) = uconvert(u"ns", inv(250u"MHz"));

    "Maximum energy DAQ can register? In keV"
    max_e::typeof(1.0*energy_unit) = 10_000u"keV" # == typemax(UInt16)

    "DAQ offset in keV"
    offset::typeof(1.0*energy_unit) = 2_000u"keV";

    "DAQ gain"
    gain = typemax(UInt16) / ((max_e+offset)/uconvert(u"keV", germanium_ionization_energy))

    trigger_threshold::Real = 0
end


function DAQmodel(sim_config::PropDict)
    if sim_config.daq.type == "generic"
        GenericDAQ(sim_config)
    else
        @info "DAQ type $(sim_config.daq.type) not implemented!"
        @info "Available types: generic"
    end
end


"""
    function GenericDAQ(sim_conf)

Construct a struct with DAQ parameters based on given simulation configuration 

sim_conf: PropDict object with simulation settings    
"""
function GenericDAQ(sim_conf::PropDict)
    GenericDAQ(
        nsamples = sim_conf.daq.nsamples,
        baseline_length = sim_conf.daq.baseline_length,
        Δt = haskey(sim_conf.daq, :sampling_rate) ? uconvert(u"ns", inv(T(sim_conf.daq.sampling_rate)u"MHz")) : T(sim_conf.daq.sampling_interval)u"ns",
        max_e = T(sim_conf.daq.max_e)u"keV",
        offset = T(sim_conf.daq.offset)u"keV",
        trigger_threshold = sim_conf.daq.trigger_threshold
    )
end

# ------------------------------------------------------------------------------------------

"""
    daq_online_filter(values, offset, window_lengths, threshold)

Simulation of the DAQ online filter (trapezoidal filter)

values: waveform y-axis values
offset: index corresponding to the first considered waveform value
window_lengths: array of 3 values corresponding to the window lengths of the trapezoidal filter
threshold: trigger threshold 

Output: float, bool
The first value is the diff returned by the trapezoidal filter.
The second value is true (triggered) or false (did not trigger)
"""
function daq_online_filter(values::AbstractVector, offset::Int, window_lengths::NTuple{3, Int}, threshold)
    wl = window_lengths
    # first window
    r1 = offset:offset+wl[1]
    # third window
    r2 = offset+wl[1]+wl[2]:offset+sum(window_lengths)
    # difference between the third and first window
    r = mean(values[r2]) - mean(values[r1])
    r, r >= threshold
end


"""
    simulate_daq(wf, daq)

Simulate DAQ effects: offset, gain, sampling rate and digitization

wf: RDWaveform object 
daq: GenericDAQ object
"""
function simulate_daq(wf::RDWaveform, daq::GenericDAQ)
    # invert the pulse if needed
    sign = integral(wf) < 0 ? -1 : 1
    wf_daq = RDWaveform(wf.time, sign * wf.value)

    # offset
    offset = uconvert(u"eV", daq.offset) / germanium_ionization_energy
    wf_daq = RDWaveform(wf_daq.time, wf_daq.value ./ustrip(germanium_ionization_energy) .+ offset)

    # gain 
    wf_daq = RDWaveform(wf_daq.time, wf_daq.value .* daq.gain)
    
    # resample -> currently simply pick every Nth sample, more elaborate simulation later
    wf_sampled = wf_daq.value[begin:Int(daq.Δt/step(wf_daq.time)):end]
    t_sampled = range(T(0)u"ns", step = daq.Δt, length = length(wf_sampled))
    wf_daq = RDWaveform(t_sampled, wf_sampled)

    # digitize
    wf_daq = RDWaveform(wf_daq.time, UInt16.(round.(wf_daq.value, digits = 0)))

    wf_daq
end



function simulate_daq(wf::RDWaveform, daq::GenericDAQ, baseline::RDWaveform)
    # invert the pulse if needed
    sign = integral(wf) < 0 ? -1 : 1
    wf_daq = RDWaveform(wf.time, sign * wf.value)

    # gain 
    wf_daq = RDWaveform(wf_daq.time, wf_daq.value .* daq.gain/ustrip(germanium_ionization_energy))
    
    # resample -> currently simply pick every Nth sample, more elaborate simulation later
    wf_sampled = wf_daq.value[begin:Int(daq.Δt/step(wf_daq.time)):end]
    t_sampled = range(T(0)u"ns", step = daq.Δt, length = length(wf_sampled))
    wf_daq = RDWaveform(t_sampled, wf_sampled)

    # extend baseline to cover all wf length
    baseline_long = extend_baseline(baseline, wf_daq)

    # now it contains offset and noise information that we need
    wf_daq = RDWaveform(wf_daq.time, wf_daq.value .+ baseline_long.value)


    # digitize
    wf_daq = RDWaveform(wf_daq.time, UInt16.(round.(wf_daq.value, digits = 0)))

    wf_daq
end


function integral(wf::RDWaveform)
    res = 0
    for i in 2:length(wf.value)
        y = (wf.value[i] - wf.value[i-1]) / 2
        dx = ustrip(wf.time[i] - wf.time[i-1])
        res += y * dx
    end
    res
end
