keV_unit = 1.0u"keV"; keV = typeof(keV_unit);
MHz_unit = 1.0u"MHz"; MHz = typeof(MHz_unit);
ns_unit = 1.0u"ns"; ns = typeof(ns_unit);
μs_unit = 1.0u"μs"; μs = typeof(μs_unit);

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

    noise_σ::Real = uconvert(u"eV", T(3)u"keV") / germanium_ionization_energy
end


"""
A shelf to put all the preamplifier parameters in
later to be filled from a configuration file
"""
@with_kw mutable struct PreAmp
    "PreAmp exp decay time"
    τ_decay::μs = T(5)*u"μs"
    
    "PreAmp rise time"
    τ_rise::ns = T(15)*u"ns"

    "sigma for electronic noise"
    noise_σ::Real = uconvert(u"eV", T(3)u"keV") / germanium_ionization_energy
end

"""
    function DAQ(elec_conf)

Construct a struct with DAQ parameters based on given electronics config PropDict
"""
function construct_DAQ(elec_conf::PropDict)
    daq = DAQ(
        nsamples=elec_conf.daq.nsamples,
        baseline_length=elec_conf.daq.baseline_length,
        sampling_rate=T(elec_conf.daq.sampling_rate)u"MHz",
        max_e=T(elec_conf.daq.max_e)u"keV",
        offset=T(elec_conf.daq.offset)u"keV",
        # sampling_rate=elec_conf.daq.sampling_rate*MHz_unit
        noise_σ=uconvert(u"eV", T(elec_conf.noise)u"keV") / germanium_ionization_energy
    )

    daq
end

"""
    function PreAmp(elec_conf)

Construct a struct with PreAmp parameters based on given electronics config PropDict
"""
function construct_PreAmp(elec_conf::PropDict)
    preamp = PreAmp(
        τ_decay=T(elec_conf.preamp.t_decay)*u"μs",
        τ_rise=T(elec_conf.preamp.t_rise)*u"ns",
        noise_σ=uconvert(u"eV", T(elec_conf.noise)u"keV") / germanium_ionization_energy
    )

    preamp
end