abstract type NoiseModel end

# add noise extracted from data
struct NoiseData <: NoiseModel
    noise_σ::Real
end

# simulate noise (fano, preamp)
struct NoiseSim <: NoiseModel
    noise_σ::Real
end

function NoiseModel(sim_config::PropDict)
    # haskey(sim_config, "noise_data") ? NoiseData(sim_config.noise_data) : NoiseSim(sim_config.preamp.noise_σ)
    # TEMPORARY
    NoiseData(sim_config.noise_data)
end



"""
    simulate_noise(wf, noise_model)

Simulate noise based on sigma extracted from data.
Should be applied after conversion to DAQ units.
Currently simple Gaussian noise.

wf: RDWaveform
noise_σ: value in DAQ units -> change to noise model 
Output: RDWaveform
"""
function simulate_noise(wf::SolidStateDetectors.RDWaveform, noise_model::NoiseModel)
    # I am no expert here. I don't know at which point one should introduce noise.
    # Also, different noise could be added at different stages. This really depends on the electronics.
    # I will just add some Gaussian Noise (σ of 3 keV defined on top)
    # lets generate 1000 random samples from this normal distribution
    gaussian_noise_dist = Normal(T(0), T(noise_model.noise_σ)) #  Normal() -> Distributions.jjl
    samples = rand(gaussian_noise_dist, 1000)

    # Now, lets add this Gaussian noise to other waveform (here, after the filters (but might be also added before))
    RDWaveform(wf.time, wf.value .+ rand!(gaussian_noise_dist, similar(wf.value)))
end