"""
Abstract type for hierarchy and multiple dispatch
that allows to choose between SolidStateDetectors and Siggen
as waveform simulation method
"""
abstract type PSSimulator end
struct SSDSimulator <: PSSimulator end
struct SiggenSimulator <: PSSimulator
    fieldgen_input::AbstractString
end


function PSSimulator(sim_config::PropDict)
    @info "Waveform simulation method: $(sim_config.sim_method)"
    if(sim_config.sim_method == "SSD")
        SSDSimulator()
    elseif(sim_config.sim_method == "siggen")
        @info "Taking fieldgen input from $(sim_config.fieldgen_input)"
        SiggenSimulator(sim_config.fieldgen_input)
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
function simulate_wf(mc_events::Table, detector_config::Dict, ::SSDSimulator)

    # simulate detector 
    det_name = detector_config["name"]
    det_h5 = joinpath("cache", basename(det_name*".h5f"))
    if isfile(det_h5)
        @info "Reading $det_name simulation from cached h5"
        simulation = SolidStateDetectors.ssd_read(det_h5, Simulation)
    else
        @info "Simulating $det_name from scratch"
        simulation = simulate_detector(detector_config)
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


function simulate_wf(mc_events::Table, detector::PropDict, siggen::SiggenSimulator)
    # to use the siggen function "read_config" from MJDSigGen.jl, we need to first create a config file it can read 

    # read lines from user input for fieldgen 
    fieldgen_lines = readlines(open(siggen.fieldgen_input))
    fieldgen_lines = replace.(fieldgen_lines, "\t" => "  ")
    # construct lines for detector related siggen input
    detector_lines = siggen_detector_config(detector)

    total_lines = vcat(["# detector related parameters"], detector_lines, fieldgen_lines)

    # write a config file 
    sig_conf_name = siggen.fieldgen_input * ".temp"
    DelimitedFiles.writedlm(sig_conf_name, total_lines)
        

    # launch siggen on our events
    setup = SigGenSetup(sig_conf_name)
    waveforms = simulate_wf(setup, mc_events)

    mcpss_table = Table(
        # channel = contact_charge_signals.chnid,
        # ievt = contact_charge_signals.evtno,
        waveform = waveforms
    )

    mcpss_mctruth = Table()

    mcpss_table, mcpss_mctruth
    
end


function simulate_wf(setup::SigGenSetup, mc_events::Table)
    # loop over events and simulate signal with siggen
    waveforms = []
    for i in 1:length(mc_events)
        push!(waveforms, simulate_wf(setup, mc_events[i].pos, mc_events[i].edep))
    end

    ArrayOfRDWaveforms(waveforms)
end


# change AbstractVector to more narrow type?
function simulate_wf(setup::SigGenSetup, pos::AbstractVector, edep::AbstractVector)
    # this is not so nice, think about it
    waveform = simulate_wf(setup, pos[1]) * ustrip(edep[1])

    for i in 2:length(pos)
        waveform = waveform .+ simulate_wf(setup, pos[i]) * ustrip(edep[i])
    end
    waveform
end

# change AbstractVector to more narrow type?
function simulate_wf(setup::SigGenSetup, pos::AbstractVector)
    # println(Tuple(round.(ustrip.(pos), digits=1)))
    # MJDSigGen.get_signal!(setup, Tuple(round.(ustrip.(pos), digits=1)))
    MJDSigGen.get_signal!(setup, Tuple(ustrip.(pos)))
end



