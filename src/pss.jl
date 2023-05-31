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

    "Path to crystal jsons"
    crystal_metadata_path::AbstractString = ""

    "Name for cached simulation. Not caching if empty string"
    cached_name::AbstractString = ""
end


"""
    SSDSimulator(sim_conf)

LegendGeSimConfig -> SSDSimulator

Construct SSDSimulator instance based on simulation
    configuration given in <sim_conf>.

Currently SSDSimulator does not have any parameters
"""
# function SSDSimulator(sim_conf::LegendGeSimConfig)
function SSDSimulator(simulation_settings::PropDict)
    coord = haskey(simulation_settings, :coordinates) ? simulation_settings.coordinates : "cylindrical"
    if !(coord in ["cartesian", "cylindrical"])
        error("$coord coordinates not implemented!\n Available: cartesian, cylindrical")
    end

    comp = haskey(simulation_settings, :computation) ? simulation_settings.computation : "2D"
    if !(comp in ["2D", "3D"])
        error("$comp computation not implemented!\n Available: 2D, 3D")
    end

    SSDSimulator(coord, comp, simulation_settings.crystal_metadata_path, simulation_settings.cached_name)
end


"""
Simulation method: siggen
"""
@with_kw struct SiggenSimulator <: PSSimulator
    "Path to fieldgen settings"
    fieldgen_config::AbstractString

    "Drift velocity correction (?)"
    drift_vel::AbstractString

    "Name for cached simulation. Not caching if empty string"
    cached_name::AbstractString = ""    
end


"""
    SiggeSimulator(sim_conf)

    LegendGeSimConfig -> SiggenSimulator

Construct SiggeSimulator instance based on simulation
    configuration given in <sim_conf>.
"""
function SiggenSimulator(simulation_settings::PropDict)
    # check if provided paths exist
    inputs = Dict(
        "fieldgen_config" => haskey(simulation_settings, :fieldgen_config) ? simulation_settings.fieldgen_config : "",
        "drift_vel" => haskey(simulation_settings, :drift_vel) ? simulation_settings.drift_vel : ""
    )
    for (param, path) in inputs
        if (path != "") && !isfile(path)
        @error "The file for $(param) that you provided does not exist: $(path)"
        end
    end

    SiggenSimulator( 
        # haskey(simulation_settings, "fieldgen_config") ? simulation_settings.fieldgen_config : "fieldgen_settings.txt",
        # haskey(simulation_settings, "drift_vel") ? simulation_settings.drift_vel : "drift_vel_tcorr.tab"
        inputs["fieldgen_config"],
        inputs["drift_vel"],
        simulation_settings.cached_name
    )
end


"""
    PSSimulator(sim_config)

LegendGeSimConfig -> <PSSimulator>

Construct a PSSSimulator supertype instance based on given simulation
    configuration <sim_config>.
Returned type depends on the simulation
    method given in the config.
"""
function PSSimulator(simulation_settings::PropDict)    
    @info "Simulation method: $(simulation_settings.method)"

    # defaults
    if(!haskey(simulation_settings, :crystal_metadata_path))
        simulation_settings[:crystal_metadata_path] = ""
        # simulation_settings.crystal_metadata_path = crystal_metadata_path
        # ToDo: move somewhere else - irrelevant if only geometry is constructed, or when simulation read from cache
        @warn "No crystal metadata path given. Simulation with dummy constant impurity density."
    elseif !ispath(simulation_settings.crystal_metadata_path)
        @error "The path to crystal metadata you provided is not valid! ($simulation_settings.crystal_metadata_path)"
    end

    if(!haskey(simulation_settings, :cached_name))
        simulation_settings[:cached_name] = ""
        @warn "No cached name was given. Not caching the simulation."
    end

    if simulation_settings.method in ["SSD", "ssd"]
        SSDSimulator(simulation_settings)
    elseif simulation_settings.method in ["siggen", "fieldgen"]
        SiggenSimulator(simulation_settings)
    else
        error("This simulation method is not implemented!")
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
            max_nsteps = 20000,
            Δt = 1u"ns",
            verbose = false);

    # SSD returns in units of "e" -> convert to eV
    n_waveforms = size(contact_charge_signals.waveform, 1)
    wf_array = Array{RDWaveform}(undef, n_waveforms)
    for i = 1:n_waveforms
        wf = contact_charge_signals.waveform[i]
        wf_array[i] = RDWaveform(wf.time, ustrip.(wf.signal) .* germanium_ionization_energy) # units eV
    end

    # ToDo: SSD returns double the number of wfs, because also inverse ones from n+ contact
    # -> just take the first half corr. to n_events in stp?
    # -> filter in a smarter way by contact?

    waveforms = ArrayOfRDWaveforms(wf_array)

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
    T = Float32 # This should be somehow defined and be passed properly
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
    waveforms = ArrayOfRDWaveforms((waveforms.time, VectorOfSimilarVectors(waveforms.signal)))    

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



