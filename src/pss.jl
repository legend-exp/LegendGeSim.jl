"""
The PSSimulator supertype defines what type of method 
    is used for the simulation of electric field and weighting potential,
    calculation of detector capacitance etc., as well as simulation of current pulses
"""
abstract type PSSimulator end


"""
Simulation method: SolidStateDetectors
"""
@with_kw struct SSDSimulator <: PSSimulator
    "Simulation coordinates (cylindrical or linear)"
    coord::AbstractString = "cylindrical"

    "Computation: 3D or 2D (phi symmetry)"
    comp::AbstractString = "2D"
end


"""
    SSDSimulator(::PropDict)

-> SSDSimulator

Construct SSDSimulator instance based on simulation
    configuration given in <sim_conf>.

Currently SSDSimulator does not have any parameters
"""
function SSDSimulator(sim_conf::PropDict)
    coord = haskey(sim_conf.simulation, :coordinates) ? sim_conf.simulation.coordinates : "cylindrical"
    if !(coord in ["cartesian", "cylindrical"])
        @error "$coord coordinates not implemented!\n Available: cartesian, cylindrical"
    end

    comp = haskey(sim_conf.simulation, :computation) ? sim_conf.simulation.computation : "2D"
    if !(comp in ["2D", "3D"])
        @error "$comp computation not implemented!\n Available: 2D, 3D"
    end

    SSDSimulator(coord, comp)
end


"""
Simulation method: siggen
"""
struct SiggenSimulator <: PSSimulator
    "Path to fieldgen settings"
    fieldgen_config::AbstractString

    "Drift velocity correction (?)"
    drift_vel::AbstractString
end


"""
    SiggeSimulator(sim_conf)

PropDict -> SiggenSimulator

Construct SiggeSimulator instance based on simulation
    configuration given in <sim_conf>.
"""
function SiggenSimulator(sim_conf::PropDict)
    @info "Taking fieldgen input from $(sim_conf.simulation.fieldgen_config)"
    SiggenSimulator(
        sim_conf.simulation.fieldgen_config,
        sim_conf.simulation.drift_vel
    )
end


"""
    PSSimulator(sim_config)

PropDict -> <PSSimulator>

Construct a PSSSimulator supertype instance based on given simulation
    configuration <sim_config>.
Returned type depends on the simulation
    method given in the config.
"""
function PSSimulator(sim_config::PropDict)
    @info "Simulation method: $(sim_config.simulation.method)"
    if sim_config.simulation.method == "SSD"
        SSDSimulator(sim_config)
    elseif sim_config.simulation.method in ["siggen", "fieldgen"]
        SiggenSimulator(sim_config)
    else
        println("This simulation method is not implemented!")
    end
end

# -------------------------------------------------------------------

"""
    simulate_waveforms(stp_events, detector)

Table, SolidStateDetectors.Simulation -> Table, Table     

Simulate pulses based on events given in <stp_events>
    using the given SSD detector simulation instance <detector>

Constructs and returns a table with resulting pulses and a pss truth table
    (may be abolished in the future as unnecessary)    
"""
function simulate_waveforms(stp_events::Table, detector::SolidStateDetectors.Simulation)
    @info("~.~.~ SolidStateDetectors")
    contact_charge_signals = SolidStateDetectors.simulate_waveforms(
            stp_events,
            detector,
            max_nsteps = 4000,
            Δt = 1u"ns",
            verbose = false);

    waveforms = ArrayOfRDWaveforms(contact_charge_signals.waveform)

    # convert to Tier1 format
    pss_table = Table(
        channel = contact_charge_signals.chnid,
        ievt = contact_charge_signals.evtno,
        waveform = waveforms
    )

    pss_truth = Table(
        detno = contact_charge_signals.detno,
        edep = contact_charge_signals.edep,
        ievt = contact_charge_signals.evtno,
        pos = contact_charge_signals.pos,
        thit = contact_charge_signals.thit
    )

    pss_table, pss_truth    
end


"""
    simulate_waveforms(stp_events, detector)

Table, AbstractString -> Table, Table     

Simulate pulses based on events given in <stp_events>
    using Siggen with settings given in <detector>

Constructs and returns a table with resulting pulses and a pss truth table
    (may be abolished in the future as unnecessary)    
"""
function simulate_waveforms(stp_events::Table, detector::SigGenSetup)
    @info("~.~.~ Siggen")

    nevents = length(stp_events)
    wf_array = Array{RDWaveform}(undef, nevents)
    for i in 1:nevents
        if(i % 100  == 0) println("$i/$nevents") end
        signal::Vector{Float32} = simulate_signal(stp_events[i].pos, stp_events[i].edep, detector)
        time = range(T(0)u"ns", step = detector.step_time_out*u"ns", length = length(signal))

        # push!(waveforms, signal)
        wf_array[i] = RDWaveform(time, signal)
    end

    ## Oliver said this should be the more efficient way of doing it
    ## as opposed to an array of separate RDWaveform instances
    # Δt = setup.step_time_out*u"ns"
    # t_end = length(waveforms[1]) * Δt
    # times = fill(0u"ns" : Δt : t_end, length(waveforms))
    # ArrayOfRDWaveforms((times, VectorOfSimilarVectors(waveforms)))

    waveforms = ArrayOfRDWaveforms(wf_array)
    # why am I doing this?
    waveforms = ArrayOfRDWaveforms((waveforms.time, VectorOfSimilarVectors(waveforms.value)))    

    pss_table = Table(
        channel = [1 for idx in 1:length(waveforms)], # lists of ADCs that triggered, 1 for HADES all the time
        ievt = stp_events.evtno,
        waveform = waveforms
    )

    # using SSD we get mc truth from SSD output
    # with siggen, let's just return the input
    # maybe mc truth will not be propagated to (sim)raw anyway
    pss_truth = Table(
        detno = stp_events.detno,
        edep = stp_events.edep,
        ievt = stp_events.evtno,
        pos = stp_events.pos,
        thit = stp_events.thit        
    )

    pss_table, pss_truth
end


"""
    simulate_wf(pos, edep, siggen_setup)

AbstractVector, AbstractVector, SigGenSetup -> Vector{Float32}

Simulate a signal from energy depositions <edep> with corresponding positions <pos>
    based on a given <siggen_setup>    
"""
function simulate_signal(pos::AbstractVector, edep::AbstractVector, siggen_setup::SigGenSetup)
    # this is not so nice, think about it
    signal = zeros(Float32, siggen_setup.ntsteps_out)

    for i in 1:length(pos)
        # in SSD the output is always in eV -> convert to eV
        signal = signal .+ MJDSigGen.get_signal!(siggen_setup, Tuple(ustrip.(pos[i]))) * ustrip(uconvert(u"eV", edep[i]))

        # somehow this doesn't work
        # MJDSigGen.get_signal!(signal, siggen_setup, Tuple(ustrip.(pos[i]))) * ustrip(uconvert(u"eV", edep[i]))
    end

    signal
end



