abstract type NoiseModel end

# add noise extracted from data
struct NoiseData <: NoiseModel
    # noise_σ::Real
    # baseline_catalog::Table
    baseline_catalog::AbstractString
    noise_σ::Real
end

# function NoiseData(baseline_file::AbstractString)
#     @info "Constructing baseline catalog for noise modeling"
#     NoiseData(baseline_catalog(baseline_file))
# end

# simulate noise (fano, preamp)
struct NoiseSim <: NoiseModel
    noise_σ::Real
end

# struct NoiseBaseline <: NoiseModel
#     baseline_database_file::AbstractString
# end

function NoiseModel(sim_config::PropDict)
    noise_model = haskey(sim_config, :noise_data) ? NoiseData(sim_config.noise_data,0) : NoiseSim(uconvert(u"eV", T(sim_config.preamp.noise_sigma)u"keV") / germanium_ionization_energy)
    # println(noise_model)
    noise_model
end


"""
    add_fnoise(mc_events, simulation, noise_model)

Add fano noise to MC events

mc_events: table in the mcstp format (output of g4_to_mcstp)
noise_model: NoiseSim object

Output: Tabless
"""
function fano_noise(mc_events::Table, det_json::AbstractString, ::NoiseSim)
    @info("Adding fano noise")
    ssd_config = detector_config(det_json, SSDSimulator())
    simulation = Simulation(SolidStateDetector{T}(ssd_config))
    det_material = simulation.detector.semiconductors[1].material
    add_fano_noise(mc_events, det_material.E_ionisation, det_material.f_fano)
    # SolidStateDetectors.add_fano_noise(mc_events, det_material.E_ionisation, det_material.f_fano)
end


function fano_noise(mc_events::Table, ::AbstractString, ::NoiseData)
    # do nothing since if we're using noise from data we do not simulate fano noise
    # not to double count
    @info("Not adding fano noise because using noise levels from data")
    mc_events
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
function simulate_noise(wf::RDWaveform, noise_model::NoiseSim)
    # I am no expert here. I don't know at which point one should introduce noise.
    # Also, different noise could be added at different stages. This really depends on the electronics.
    # I will just add some Gaussian Noise (σ of 3 keV defined on top)
    # lets generate 1000 random samples from this normal distribution
    gaussian_noise_dist = Normal(T(0), T(noise_model.noise_σ)) #  Normal() -> Distributions.jjl

    # Now, lets add this Gaussian noise to other waveform (here, after the filters (but might be also added before))
    RDWaveform(wf.time, wf.value .+ rand!(gaussian_noise_dist, similar(wf.value)))
end

