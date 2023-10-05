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
    # Note: these default T() expressions don't work cause T not defined
    "PreAmp exp decay time"
    τ_decay::typeof(1.0*μs_unit) = T(5)*u"μs"
    
    "PreAmp rise time"
    τ_rise::typeof(1.0*ns_unit) = T(15)*u"ns"

    "Preamp offset in keV"
    offset_keV::typeof(1.0*energy_unit) = 0u"keV";

    "Preamp offset in ADC (alternative)"
    offset_ADC = 0; # ADC units

    # "Maximum DAQ energy"
    # max_e::typeof(1.0*energy_unit) = 0u"keV" # == typemax(UInt16)

    "Preamp gain"
    # If offset = 0, it means gain has to be inferred (from the data baseline offset)
    # How does it end up being Float64?
    # gain =             
    #     if offset_keV == 0u"keV"
    #         if offset_ADC == 0
    #             0
    #         else
    #             # ( typemax(UInt16) - offset_ADC ) * uconvert(u"keV", germanium_ionization_energy) / max_e # ADC/charge
    #             ( typemax(UInt16) - offset_ADC ) / uconvert(u"eV", max_e) # ADC/eV
    #         end
    #     else
    #         # typemax(UInt16) / ((max_e+offset_keV)/uconvert(u"keV", germanium_ionization_energy)) # ADC/charge
    #         typemax(UInt16) / uconvert(u"eV", max_e + offset_keV) # ADC/eV
    #     end            
    gain::typeof(1.0/1u"eV")
    # gain = offset == 0u"keV" ? 0 : typemax(UInt16) / ((max_e+offset)/uconvert(u"keV", germanium_ionization_energy))
    # gain = offset_ADC == 0 ? : ( typemax(UInt16) - offset_ADC ) * uconvert(u"keV", germanium_ionization_energy) / max_e

    "Preamp Gaussian noise in ADC"
    noise_σ_ADC = 0

    "Preamp Gaussian noise in keV"
    noise_σ_keV::typeof(1.0*energy_unit)= 0u"keV"
    # noise_σ_keV::typeof(1.0*energy_unit) = noise_σ_ADC / gain

end


"""
CC4 circuit model (ToDo)
"""
struct CC4 <: PreAmp end


"""
    GenericPreAmp(sim_conf)

LegendGeSimConfig -> GenericPreAmp    

Construct a GenericPreAmp instance based on given simulation configuration <sim_conf>
"""
function GenericPreAmp(preamp_config::PropDict)
    T = Float32 # This should be somehow defined and be passed properly
    GenericPreAmp(
        τ_decay=T(preamp_config.t_decay)*u"μs",
        τ_rise=T(preamp_config.t_rise)*u"ns",
        # optional, if not given usually means NoiseFromData is used -> make a better system for managing this!
        offset_keV = haskey(preamp_config, :offset_in_keV) ? T(preamp_config.offset_in_keV)u"keV" : 0u"keV",
        offset_ADC = haskey(preamp_config, :offset_in_ADC) ? T(preamp_config.offset_in_ADC) : 0,
        # max_e = T(preamp_config.max_e)u"keV",
        gain = preamp_config.gain/u"eV",
        # optional, if not given usually means NoiseFromData is used -> make a better system for managing this!
        noise_σ_ADC = haskey(preamp_config, :noise_sigma_ADC) ? T(preamp_config.noise_sigma_ADC) : 0,
        noise_σ_keV = haskey(preamp_config, :noise_sigma_keV) ? T(preamp_config.noise_sigma_keV)u"keV" : 0u"keV"
    )
    # ToDo ! insert check if both offset in keV and ADC is provided is bad
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
        GenericPreAmp(preamp_settings)
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

    ## noise before preamp
    # Alessandro Razeto said noise should be added before preamp
    # if I simply move adding Gaus noise here, it gives weird cyclical vibrations in the waveform
    # then after build_dsp there is no Amax at all, I guess DSP fails
    # wf_noise = simulate_noise(wf, preamp)

    ## rise and decay time
    csa_filter = RadiationDetectorDSP.simple_csa_response_filter(
        preamp.τ_rise / step(wf.time),
        uconvert(u"ns", preamp.τ_decay) / step(wf.time))

    wf_preamp = RDWaveform(wf.time, filt(csa_filter, wf.signal)) 

    wf_preamp = simulate_noise(wf_preamp, preamp)

    ## offset
    # wf values are in eV but without u"eV" attached
    offset_keV = ustrip(uconvert(u"eV", preamp.offset_keV))
    # if ADC offset is used, offset_keV is zero
    wf_preamp = RDWaveform(wf_preamp.time, wf_preamp.signal .+ offset_keV)

    ## gain 
    # gain is in eV^-1
    gain = ustrip(preamp.gain)
    wf_preamp = RDWaveform(wf_preamp.time, wf_preamp.signal * gain)

    ## ADC offset: if keV offset is used, ADC offset is zero
    wf_preamp = RDWaveform(wf_preamp.time, wf_preamp.signal .+ preamp.offset_ADC)
      
    ## what is this step? -> now gain is in ADC/eV, current was in eV, final result is in ADC
    # wf_preamp = RDWaveform(wf_preamp.time, wf_preamp.signal ./ ustrip(germanium_ionization_energy))

    wf_preamp
end