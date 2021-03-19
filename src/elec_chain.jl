# keV_unit = 1.0u"keV"; keV = typeof(keV_unit);
# MHz_unit = 1.0u"MHz"; MHz = typeof(MHz_unit);
# const freq_unit = u"MHz"
# ns_unit = 1.0u"ns"; ns = typeof(ns_unit);
# μs_unit = 1.0u"μs"; μs = typeof(μs_unit);

ns_unit = u"ns"
μs_unit = u"μs"

T = Float32
germanium_ionization_energy = T(2.95)u"eV"

abstract type ElecChain end

"""
A shelf to put all the preamplifier parameters in
later to be filled from a configuration file
"""
@with_kw struct PreAmp <: ElecChain
    "PreAmp exp decay time"
    τ_decay::typeof(1.0*μs_unit) = T(5)*u"μs"
    
    "PreAmp rise time"
    τ_rise::typeof(1.0*ns_unit) = T(15)*u"ns"

    # "sigma for electronic noise"
    # noise_σ::Real = uconvert(u"eV", T(3)u"keV") / germanium_ionization_energy
end


"""
    function PreAmp(elec_conf)

Construct a struct with PreAmp parameters based on given electronics config PropDict
"""
function construct_PreAmp(sim_conf::PropDict)
    PreAmp(
        τ_decay=T(sim_conf.preamp.t_decay)*u"μs",
        τ_rise=T(sim_conf.preamp.t_rise)*u"ns"
        # noise_σ=uconvert(u"eV", T(elec_conf.noise)u"keV") / germanium_ionization_energy
    )
end


"""
    simulate_elec(wf, preamp)

Simulate a charge sensitive amplifier (`CSA`)
Here, the parameters τ_rise and τ_decay have to be given in units of samples,
because the `filt`-function does not know the Δt between the samples.

wf: RDWaveform
preamp: a PreAmp struct
Output: RDWaveform
"""
function simulate_elec(wf::SolidStateDetectors.RDWaveform, preamp::PreAmp)
    csa_filter = dspjl_simple_csa_response_filter(
        preamp.τ_rise / step(wf.time),
        uconvert(u"ns", preamp.τ_decay) / step(wf.time))

    RDWaveform(wf.time, filt(csa_filter, wf.value))
end