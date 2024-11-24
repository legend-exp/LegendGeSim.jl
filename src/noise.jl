"""
The NoiseModel supertype corresponds to different methods
    of simulating noise coming from the electronics chain
    components.

It is not simply a string in the ElecChain struct because
    it must exist as a separate instance for the stp->pss step 
    that is independent from the elec simulation, but has to know
    what type of noise model is used to deal with fano noise.
"""
abstract type NoiseModel end

# in case no noise wanted just ideal pulses
struct NoiseNone <: NoiseModel end

"""
The NoiseFromSim model means simulating noise from scratch
    starting from fano noise of the germanium crystal
    and ending with noise coming from electronics components.
"""
struct NoiseFromSim <: NoiseModel end
# @with_kw struct NoiseFromSim <: NoiseModel
#     "Noise σ in keV for Gaussian noise simulation"
#     noise_σ::typeof(1.0*energy_unit) = 0u"keV"
# end



"""
The NoiseFromData model means instead of simulating given 
    preamp noise, and also gain and offset, these parameters will be
    accounted for or inferred from real data baselines
    contained in the <baseline_catalog>
"""
struct NoiseFromData <: NoiseModel
    "A table of data baselines in RDWaveform format"
    baseline_catalog::Table
end


"""
    NoiseFromData(sim_conf)

LegendGeSimConfig -> NoiseFromData 

Construct a NoiseFromData struct based on given simulation configuration <sim_conf>
"""
function NoiseFromData(noise_settings::PropDict)
    @info "//\\//\\// Noise levels and offset added via slapping the baseline from data on top of the waveform"
    # construct table of baselines based on given raw data hdf5 file
    baseline_table = baseline_catalog(noise_settings.path_to_data_file)
    NoiseFromData(baseline_table) 
end


"""
    NoiseModel(sim_config)

LegendGeSimConfig -> <NoiseModel>

Constuct a NoiseModel supertype instance based on simulation settings 
    given in <sim_config>
Type of <NoiseModel> depends on <sim_config> settings.
"""
function NoiseModel(noise_settings::PropDict)
    # if !haskey(noise_settings, :type)
    if noise_settings.type == "none"
        NoiseNone()
    elseif noise_settings.type == "sim"
        NoiseFromSim()
    else
        NoiseFromData(noise_settings)
    end
    # ToDo add elseif and error if not none, sim or data
end

# -------------------------------------------------------------------

"""
    fano_noise(events, det_meta, env, ::NoiseFromSim)

Table, PropDict, Environment -> Table 

Calculate fano noise level based on the detector specification provided in
    LEGEND metadata <det_meta> and environment settings provided in <env>,
    and add it to given <events>
"""
function fano_noise(events::Table, det_meta::PropDict, env::Environment, ::NoiseFromSim)
    @info "//\\//\\//\\ Adding fano noise"
    detector = LEGEND_SolidStateDetector(Float32, det_meta, env)
    det_material = detector.semiconductor.material
    add_fano_noise(events, det_material.E_ionisation, det_material.f_fano)
end


"""
    fano_noise(events, ::PropDict, ::Environment, ::NoiseFromData)

Table -> Table    

Do nothing since we do not need to simulate fano noise separately when 
    using data baselines to account for noise levels.
"""
function fano_noise(events::Table, ::PropDict, ::Environment, ::NoiseFromData)
    println("Not adding fano noise because using noise levels from data")
    events
end

function fano_noise(events::Table, ::PropDict, ::Environment, ::NoiseNone)
    println("No fano noise (no noise model given)")
    events
end

"""
    simulate_noise(wf, preamp)

RDWaveform, PreAmp -> RDWaveform

Simulate effects of the preamplifier <preamp> on the given waveform <wf>.
"""
function simulate_noise(wf::RDWaveform, preamp::PreAmp)
    T = Float32 # This should be somehow defined and be passed properly
    # if we're not simulating noise from scratch, we'll do this anyway, but noise sigma will be 0
    # -> maybe there's a better way?
    # -> also will double count if user puts NoiseNone or NoiseFromData but forgets a non-zero sigma

    # wf values are in eV (without u"eV" units attached), noise sigma is in keV
    noise_σ::T = T(ustrip(u"eV", iszero(preamp.noise_σ_keV) ? preamp.noise_σ_ADC / preamp.gain : uconvert(u"eV", preamp.noise_σ_keV)))
    
    # --- white Gaus
    # gaussian_noise_dist = Normal(T(0), T(noise_σ))
    # RDWaveform(wf.time, wf.signal .+ rand!(gaussian_noise_dist, similar(wf.signal)))
    # --- pink Gaus
    gaus_pink = SignalAnalysis.PinkGaussian(length(wf.signal), noise_σ)
    RDWaveform(wf.time, wf.signal .+ rand(gaus_pink))
end


"""
    simulate_noise(wf, noise_model)

RDWaveform, NoiseFromData -> RDWaveform

Simulate noise and offset by picking a random baseline from the baseline catalog
    contained in <noise_model> and slapping it on top of the given waveform <wf>.
"""
function simulate_noise(wf::RDWaveform, noise_model::NoiseFromData)
    # NOTE: has to be done after FADC i.e. after ElecChain
    baseline = rand(Tables.getcolumn(Tables.columns(noise_model.baseline_catalog), :waveform))
    # extend to match the wf
    baseline_long = extend_baseline(baseline, wf)
    # slap on top of the waveform as promised
    RDWaveform(wf.time, wf.signal .+ baseline_long.signal)
end


