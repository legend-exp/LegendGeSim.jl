"""
Abstract type for hierarchy and multiple dispatch
that allows to choose between SolidStateDetectors and Siggen
as waveform simulation method
"""
abstract type PSSimulator end
@with_kw struct SSDSimulator <: PSSimulator
    crystal_t::Real = 0
end

struct SiggenSimulator <: PSSimulator
    crystal_t::Real
    fieldgen_input::AbstractString
end



function PSSimulator(sim_config::PropDict)
    @info "Waveform simulation method: $(sim_config.sim_method)"
    if(sim_config.sim_method == "SSD")
        SSDSimulator(sim_config.crystal_t)
    elseif(sim_config.sim_method == "siggen")
        @info "Taking fieldgen input from $(sim_config.fieldgen_input)"
        SiggenSimulator(sim_config.crystal_t, sim_config.fieldgen_input)
    else
        println("This simulation method is not implemented!")
    end
end



"""
    simulate_wf(mc_events, simulation, simulator)

Simulate waveforms based on given MC events and detector simulation,
contained in the MCstp struct. Returns resulting mcpss table
and a secon table containing MC truth

mcstp: MCstp struct
simulator: Simulator object inheriting from PSSimulator for multiple dispatch

Output: Table, Table
"""
function simulate_wf(mc_events::Table, det_json::AbstractString, ssd_sim::SSDSimulator)

    # simulate detector 
    det_name = splitext(basename(det_json))[1]
    det_h5 = joinpath("cache", det_name*".h5f")
    if isfile(det_h5)
        @info "Reading $det_name simulation from cached h5"
        simulation = SolidStateDetectors.ssd_read(det_h5, Simulation)
    else
        @info "Simulating $det_name from scratch"
        ssd_config = detector_config(det_json, ssd_sim)
        simulation = simulate_detector(ssd_config)
    end

    @info("Simulating waveforms")
    contact_charge_signals = SolidStateDetectors.simulate_waveforms(
            mc_events,
            simulation,
            max_nsteps = 4000,
            Δt = 1u"ns",
            verbose = false);

    waveforms = ArrayOfRDWaveforms(contact_charge_signals.waveform)

    # convert to Tier1 format
    mcpss_table = Table(
        channel = contact_charge_signals.chnid,
        ievt = contact_charge_signals.evtno,
        waveform = waveforms
    )

    mcpss_mctruth = Table(
        detno = contact_charge_signals.detno,
        edep = contact_charge_signals.edep,
        ievt = contact_charge_signals.evtno,
        pos = contact_charge_signals.pos,
        thit = contact_charge_signals.thit
    )

    mcpss_table, mcpss_mctruth
end


function simulate_wf(mc_events::Table, det_json::AbstractString, siggen_sim::SiggenSimulator)
    # to use the siggen function "read_config" from MJDSigGen.jl, we need to first create a config file it can read 

    siggen_setup = detector_config(det_json, siggen_sim)
    waveforms = simulate_wf(siggen_setup, mc_events)

    mcpss_table = Table(
        channel = [1 for idx in 1:length(waveforms)], # lists of ADCs that triggered, 1 for HADES all the time
        ievt = mc_events.evtno,
        waveform = waveforms
    )

    # using SSD we get mc truth from SSD output
    # with siggen, let's just return the input
    # maybe mc truth will not be propagated to mcraw anyway
    mcpss_mctruth = Table(
        detno = mc_events.detno,
        edep = mc_events.edep,
        ievt = mc_events.evtno,
        pos = mc_events.pos,
        thit = mc_events.thit        
    )

    mcpss_table, mcpss_mctruth
    
    
end


function simulate_wf(setup::SigGenSetup, mc_events::Table)
    # loop over events and simulate signal with siggen
    waveforms = Vector{Vector{Float32}}()
    nevents = length(mc_events)
    for i in 1:nevents
        if(i % 100  == 0) println("$i/$nevents") end
        signal::Vector{Float32} = simulate_wf(setup, mc_events[i].pos, mc_events[i].edep)
        push!(waveforms, signal)
    end

    Δt = setup.step_time_out*u"ns"
    t_end = length(waveforms[1]) * Δt
    times = fill(0u"ns" : Δt : t_end, length(waveforms))

    ArrayOfRDWaveforms((times, VectorOfVectors(waveforms)))

end


# change AbstractVector to more narrow type?
function simulate_wf(setup::SigGenSetup, pos::AbstractVector, edep::AbstractVector)
    # this is not so nice, think about it
    signal = simulate_wf(setup, pos[1]) * ustrip(edep[1])

    for i in 2:length(pos)
        signal = signal .+ simulate_wf(setup, pos[i]) * ustrip(edep[i])*1e6 # TEMP
    end
    signal
end

# change AbstractVector to more narrow type?
function simulate_wf(setup::SigGenSetup, pos::AbstractVector)
    MJDSigGen.get_signal!(setup, Tuple(ustrip.(pos)))
end



