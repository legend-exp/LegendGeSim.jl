"""
The ElecChain supertype corresponds to the chain of 
    electronic components involved in the DAQ setup
    that appear before the trigger.
"""
abstract type ElecChain end

"""
Generic electronics chain consisting of a preamplifier
    and an FADC module
"""
mutable struct GenericElecChain <: ElecChain
    preamp::PreAmp

    fadc::FADC
end


"""
    GenericElecChain(sim_conf)

LegendGeSimConfig -> GenericElecChain

Construct electronics components based on simulation 
    configuration <sim_conf> and create a GenericChain instance
    based on these components.
"""
function GenericElecChain(setup_settings::PropDict)
    preamp = PreAmp(setup_settings.preamp)
    fadc = FADC(setup_settings.fadc)

    GenericElecChain(preamp, fadc)
end


"""
    ElecChain(sim_conf)

LegendGeSimConfig -> <ElecChain>

Construct an ElecChain supertype struct based on given simulation configuration.
Type of returned instance depends on settings in <sim_conf>
Currently only one type of ElecChain available (GenericElecChain),
    rendering this function temporarily redundant.
"""
function ElecChain(setup_settings::PropDict)
    GenericElecChain(setup_settings)
end


# -------------------------------------------------------------------

"""
    simulate(wf, elec_chain, ::NoiseFromSim)

RDWaveform, ElecChain -> RDWaveform

Simulate effects of the electronics chain on the waveform
    modelling noise based on each component.

"""
function simulate(wf::RDWaveform, elec_chain::ElecChain, ::NoiseFromSim)
    # simply simulate elec chain the way it is
    # noise simulation of components is included in the simulation of components
    simulate(wf, elec_chain)
end


"""
    simulate(wf, elec_chain, noise_model)

RDWaveform, GenericElecChain, NoiseFromData -> RDWaveform

Simulate effects of the electronics chain on the waveform
    modelling noise based on baselines extracted from data
    (note: also takes preamp offset into account)
"""
function simulate(wf::RDWaveform, elec_chain::GenericElecChain, noise_model::NoiseFromData)
    ## simulate elec chain

    # check if noise and offset are double counted
    if (elec_chain.preamp.noise_σ != 0u"keV") || (elec_chain.preamp.offset != 0u"keV")
        @info "WARNING!\n
        You are simulating noise and offset based on data baselines,
        but your PreAmp settings are:\n
        -- noise sigma = $(elec_chain.preamp.noise_σ)\n
        -- offset = $(elec_chain.preamp.offset)"
        @info "Proceeding with the simulation, but there might be double counting"
    end

    # simulate elec chain
    wf_elec = simulate(wf, elec_chain)

    ## simulate noise and offset via slapping data baselines on top of the resulting waveform
    simulate_noise(wf_elec, noise_model)
end


"""
    simulate(wf, elec_chain)

RDWaveform, ElecChain -> RDWaveform

Simulate effects of the electronics chain components on the waveform.
"""
function simulate(wf::RDWaveform, elec_chain::ElecChain)
    # elec chain components
    component_names = fieldnames(typeof(elec_chain))
    
    # start with originally given waveform
    wf_elec = wf
    # simulate effect of each component
    for name in component_names
        component = getfield(elec_chain, name)
        wf_elec = simulate(wf_elec, component)
    end

    wf_elec
end

    