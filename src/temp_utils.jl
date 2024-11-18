# this should exist?
"""
    filename(path)

Sring -> String

Extract core name of the file from path 

E.g. filename("/some/path/to/file_name.ext") -> "file_name"
"""
filename(path) = splitext(basename(path))[1]


"""
    remove_negative(value)

Real -> Real 

Replace negative values by zeros.

Used to remove a glitch in SSD pulses as a quick fix
    while the glitch is being debugged.
"""
function remove_negative(value::T)::T where {T <: Real}
    max(zero(T), value)
end


"""
    preamp_gain(preamp, noise_model)

GenericPreAmp, NoiseFromData -> Real

Calculate gain of <preamp> based on its max energy and 
    the offset observed in data baselines contained in <noise_model>

In principle we should not do this in this simulation, and user has to provide 
    the precise parameters of the electronics chain used in producing data 
    the baselines from which are contained in <noise_model>.
"""
function preamp_gain(preamp::PreAmp, noise_model::NoiseFromData)
    # average offset in baselines 
    offset = mean(mean.(noise_model.baseline_catalog.waveform.signal))
    # calculate gain based on offset in data baselines 
    (typemax(UInt16) - offset) * germanium_ionization_energy / uconvert(u"eV", preamp.max_e)
end


"""
    preamp_gain(preamp, ::NoiseFromSim)

GenericPreAmp -> Float64

Do nothing and return original value of gain parameter in <preamp>,
    since when NoiseFromSim model is used, gain is provided by the user
    (or calculated based on user given max energy)    
"""
function preamp_gain(preamp::GenericPreAmp, ::Union{NoiseFromSim,NoiseNone})
    preamp.gain
end


# """
#     trigger_threshold(trigger, noise_model, preamp)

# Trigger, NoiseFromSim, GenericPreAmp -> Real 

# Calculate trigger threshold based on whether
#     the user provided it in the simulation settings
# """
# function trigger_threshold(trigger::Trigger, noise_model::NoiseFromSim, preamp::PreAmp)
#     # if trigger threshold not given, take noise level as reference 
#     threshold = trigger.threshold_keV == 0u"keV" ? noise_model.noise_σ * 3 : trigger.threshold_keV
#     uconvert(u"eV", threshold) / germanium_ionization_energy * preamp.gain
# end 


"""
    trigger_threshold(trigger, preamp, ::NoiseFromSim)

Trigger, GenericPreAmp -> Real 

In NoiseFromSim setting, non-zero trigger threshold in keV MUST be given by the user.
Calculate final threshold based the value contained in <trigger>, and <preamp> gain.

"""
function trigger_threshold(trigger::Trigger, preamp::GenericPreAmp, ::Union{NoiseFromSim, NoiseNone})
    # uconvert(u"eV", trigger.threshold_keV) / germanium_ionization_energy * preamp.gain
    uconvert(u"eV", trigger.threshold_keV) * preamp.gain # gain now in ADC/eV
    # threshold = trigger.threshold_keV == 0u"keV" ? noise_model.noise_σ * 3 : trigger.threshold_keV
    # uconvert(u"eV", threshold) / germanium_ionization_energy * preamp.gain
end 


"""
    trigger_threshold(trigger, preamp, noise_model)

Trigger, GenericPreAmp, NoiseFromData -> Real 

In NoiseFromData setting, if trigger threshold in keV is not provided,
    it is calculated based on noise levels in the baselines contained in
    <noise_model>.
Otherwise the final threshold is calculated based on <preamp> gain.
"""
function trigger_threshold(trigger::Trigger, preamp::GenericPreAmp, noise_model::NoiseFromData)
    trigger.threshold_keV == 0u"keV" ? std(noise_model.baseline_catalog.waveform[1].signal) * 3 : uconvert(u"eV", trigger.threshold_keV) / germanium_ionization_energy * preamp.gain
end

capacitance_matrix(sim::Simulation) = calculate_capacitance_matrix(sim)