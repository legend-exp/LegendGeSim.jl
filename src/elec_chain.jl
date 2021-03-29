ns_unit = u"ns"
μs_unit = u"μs"
T = Float32


# germanium_ionization_energy = T(2.95)u"eV"
germanium_ionization_energy = SolidStateDetectors.material_properties[:HPGe].E_ionisation # already in eV


"""
Abstract Electronics Chain type for hierarchy and multiple dispatch
"""
abstract type ElecChain end


"""
A shelf to put all the preamplifier parameters in
"""
@with_kw struct PreAmp <: ElecChain
    "PreAmp exp decay time"
    τ_decay::typeof(1.0*μs_unit) = T(5)*u"μs"
    
    "PreAmp rise time"
    τ_rise::typeof(1.0*ns_unit) = T(15)*u"ns"

    "sigma for electronic noise"
    noise_σ::Real = uconvert(u"eV", T(3)u"keV") / germanium_ionization_energy
end


"""
    function PreAmp(sim_conf_file)

Construct a struct with Preamp parameters based on given simulation configuration file

sim_conf_file: json file with simulation settings    
"""
function PreAmp(sim_conf::AbstractString)
    PreAmp(LegendGeSim.load_config(sim_conf))
end


"""
    function PreAmp(sim_conf)

Construct a struct with PreAmp parameters based on given simulation configuration 

sim_conf: PropDict object with simulation settings    
"""
function PreAmp(sim_conf::PropDict)
    PreAmp(
        τ_decay=T(sim_conf.preamp.t_decay)*u"μs",
        τ_rise=T(sim_conf.preamp.t_rise)*u"ns",
        noise_σ = haskey(sim_conf.preamp, :noise_sigma) ? uconvert(u"eV", T(sim_conf.preamp.noise_sigma)u"keV") / germanium_ionization_energy : 0
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
function simulate_elec(wf::RDWaveform, preamp::PreAmp)
    csa_filter = dspjl_simple_csa_response_filter(
        preamp.τ_rise / step(wf.time),
        uconvert(u"ns", preamp.τ_decay) / step(wf.time))

    RDWaveform(wf.time, filt(csa_filter, wf.value))
end