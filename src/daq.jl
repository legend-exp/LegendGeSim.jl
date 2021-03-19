# keV_unit = 1.0u"keV"; keV = typeof(keV_unit);
const energy_unit = u"keV"

"""
Abstract DAQ type for hierarchy
"""
abstract type DAQ
end

"""
A shelf to put all the DAQ parameters in
later to be filled from a configuration file
"""
@with_kw struct GenericDAQ <: DAQ
    "How many samples DAQ stores (final wf length). Has to be smaller than total wf length (defined in mcraw_to_mcpss.jl)"
    nsamples::Int = 4000; # samples per waveform

    "Length of DAQ baseline in number of samples"
    baseline_length::Int = 400;

    # "DAQ sampling rate. Needed to calculate Δt"
    # sampling_rate::MHz = 250u"MHz"

    "Inverse of sampling rate. Needed for DAQ trigger"
    Δt::typeof(1.0*ns_unit) = uconvert(u"ns", inv(250u"MHz"));

    # "DAQ type. Needed to make sure the wf values are..."
    # daq_type = UInt16

    "Maximum energy DAQ can register? In keV"
    # max_e::keV = 10_000u"keV" # == typemax(UInt16)
    max_e::typeof(1.0*energy_unit) = 10_000u"keV" # == typemax(UInt16)

    "DAQ offset in keV"
    # offset::keV = 2_000u"keV";
    offset::typeof(1.0*energy_unit) = 2_000u"keV";

    "What is this parameter? -> DAQ gain, no units"
    gain = typemax(UInt16) / ((max_e+offset)/uconvert(u"keV", germanium_ionization_energy))

    # noise in DAQ units?.. needed to compute trigger threshold
    noise_σ = uconvert(u"eV", T(3)u"keV") / germanium_ionization_energy
    
    trigger_threshold::Real = noise_σ * 10 * gain
end


"""
    function DAQ(elec_conf)

Construct a struct with DAQ parameters based on given electronics config PropDict
"""
function construct_GenericDAQ(sim_conf::PropDict)
    GenericDAQ(
        nsamples=sim_conf.daq.nsamples,
        baseline_length=sim_conf.daq.baseline_length,
        # sampling_rate=T(elec_conf.daq.sampling_rate)u"MHz",
        Δt=uconvert(u"ns", inv(T(sim_conf.daq.sampling_rate)u"MHz")),
        max_e=T(sim_conf.daq.max_e)u"keV",
        offset=T(sim_conf.daq.offset)u"keV",
        # sampling_rate=elec_conf.daq.sampling_rate*MHz_unit
        # noise_σ=uconvert(u"eV", T(sim_conf.daq.noise)u"keV") / germanium_ionization_energy
        # daq_trigger_threshold = sim_conf.noise_data * 10 * daq.c
        noise_σ = sim_conf.noise_data
    )
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


function simulate_daq(wf::SolidStateDetectors.RDWaveform, daq::GenericDAQ)
    # offset
    o = daq.gain * uconvert(u"eV", daq.offset) / germanium_ionization_energy

    # invert the pulse if needed
    sign = wf.value[end] < 0 ? -1 : 1

    # apply gain and offset
    # daq_buffer_wv = RDWaveform(wf.time, daq.daq_type.(round.(sign * wf.value .* daq.gain .+ o, digits = 0)))
    # daq_buffer_wv = RDWaveform(wf.time, UInt16.(round.(sign / ustrip(germanium_ionization_energy) * wf.value .* daq.gain .+ o, digits = 0)))
    wf_off = RDWaveform(wf.time, sign / ustrip(germanium_ionization_energy) * wf.value .* daq.gain .+ o)

    ts = range(T(0)u"ns", step = daq.Δt, length = daq.nsamples)
    # resample -> currently simply pick every Nth sample, more elaborate simulation later
    wf_samp = RDWaveform(ts, wf_off.value[range(1, step = Int(daq.Δt/step(wf.time)), length = daq.nsamples)])
    # digitize
    wf_dig = RDWaveform(wf_samp.time, UInt16.(round.(wf_samp.value, digits = 0)))

    wf_dig
end