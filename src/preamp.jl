"""
The PreAmp supertype corresponds to the charge sensitive preamplifier component
"""
abstract type PreAmp end


"""
The GenericPreAmp is a dummy PreAmp model that accounts for
    effects such as decay and rise time, offset, noise and gain.

This dummy model involves no real electronics response, only gain.

This struct is currently mutable because in the case of noise modelling based on data,
    the gain has to match the one in the baselines extracted from data,
    so in that case the PreAmp is initialized with gain = 0, and the gain
    is calculated later. I assume this is temporary, as the user should know which
    gain / max_e the data was produced with and give it as simulation input.
"""
@with_kw mutable struct GenericPreAmp <: PreAmp
    "PreAmp exp decay time"
    τ_decay::typeof(1.0*μs_unit) = T(5)*u"μs"
    
    "PreAmp rise time"
    τ_rise::typeof(1.0*ns_unit) = T(15)*u"ns"

    "Preamp offset in keV"
    offset::typeof(1.0*energy_unit) = 0u"keV";

    "Maximum DAQ energy"
    max_e::typeof(1.0*energy_unit) = 10_000u"keV" # == typemax(UInt16)

    "Preamp gain"
    # If offset = 0, it means gain has to be inferred (from the data baseline offset)
    # How does it end up being Float64?
    gain = offset == 0u"keV" ? 0 : typemax(UInt16) / ((max_e+offset)/uconvert(u"keV", germanium_ionization_energy))

    "Preamp Gaussian noise"
    noise_σ::typeof(1.0*energy_unit)= 0u"keV"
end


"""
CC2 circuit model (ToDo)
what. it's CC4. what is this
"""
struct CC2 <: PreAmp end


"""
    GenericPreAmp(sim_conf)

LegendGeSimConfig -> GenericPreAmp    

Construct a GenericPreAmp instance based on given simulation configuration <sim_conf>
"""
# function GenericPreAmp(sim_conf::LegendGeSimConfig)
#     T = Float32 # This should be somehow defined and be passed properly
#     GenericPreAmp(
#         τ_decay=T(sim_conf.dict.setup.preamp.t_decay)*u"μs",
#         τ_rise=T(sim_conf.dict.setup.preamp.t_rise)*u"ns",
#         offset = T(sim_conf.dict.setup.preamp.offset)u"keV",
#         max_e = T(sim_conf.dict.setup.preamp.max_e)u"keV",
#         noise_σ = haskey(sim_conf.dict.setup.preamp, :noise_sigma) ? T(sim_conf.dict.setup.preamp.noise_sigma)u"keV" : 0u"keV"
#     )
# end

function GenericPreAmp(preamp_config::PropDict)
    T = Float32 # This should be somehow defined and be passed properly
    GenericPreAmp(
        τ_decay=T(preamp_config.t_decay)*u"μs",
        τ_rise=T(preamp_config.t_rise)*u"ns",
        offset = T(preamp_config.offset)u"keV",
        max_e = T(preamp_config.max_e)u"keV",
        noise_σ = haskey(preamp_config, :noise_sigma) ? T(preamp_config.noise_sigma)u"keV" : 0u"keV"
    )
end


"""
    PreAmp(sim_conf)

LegendGeSimConfig -> <PreAmp>    

Construct a PreAmp supertype instance based on given simulation configuration <sim_conf>.
Returned type depends on the settings in <sim_conf>.
Currently only GenericPreAmp is available.
"""
function PreAmp(preamp_settings::PropDict)
    if preamp_settings.type == "generic"
        GenericPreamp(preamp_settings)
    else
        error("Preamp type $(preamp_settings.type) not implemented!\n
        Available type: generic")
    end
end

# function PreAmp(sim_conf::LegendGeSimConfig)
#     GenericPreAmp(sim_conf)
# end

# ------------------------------------------------------------------------------------------

"""
    simulate(wf, preamp)

RDWaveform, GenericPreAmp -> RDWaveform

Simulate the effecs of the preamp on the waveform.
"""
function simulate(wf::RDWaveform, preamp::GenericPreAmp)
    ## rise and decay time
    csa_filter = RadiationDetectorDSP.simple_csa_response_filter(
        preamp.τ_rise / step(wf.time),
        uconvert(u"ns", preamp.τ_decay) / step(wf.time))

    wf_preamp = RDWaveform(wf.time, filt(csa_filter, wf.signal))

    ## noise
    wf_preamp = simulate_noise(wf_preamp, preamp) 

    ## offset
    # wf values are in eV but without u"eV" attached
    offset = ustrip(uconvert(u"eV", preamp.offset))
    wf_preamp = RDWaveform(wf_preamp.time, wf_preamp.signal .+ offset)

    ## gain 
    wf_preamp = RDWaveform(wf_preamp.time, wf_preamp.signal * preamp.gain)
      
    ## what is this step?
    # normalizing by germanium ionization energy?
    # basically converting to number of charge carriers?
    # but after the gain?
    # I guess it's like a dummy MeV -> ADC (calibration curve?)
    wf_preamp = RDWaveform(wf_preamp.time, wf_preamp.signal ./ ustrip(germanium_ionization_energy))

    wf_preamp
end