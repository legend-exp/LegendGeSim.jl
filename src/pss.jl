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

    time_step::typeof(1.0*ns_unit) = 1u"ns"

    diffusion::Bool = false

    self_repulsion::Bool = false

    number_of_carriers::Int = 1
end


"""
    SSDSimulator(sim_conf)

LegendGeSimConfig -> SSDSimulator

Construct SSDSimulator instance based on simulation
    configuration given in <sim_conf>.

Currently SSDSimulator does not have any parameters
"""
@with_kw struct SiggenSimulator <: PSSimulator
    "Path to fieldgen settings"
    fieldgen_config::AbstractString

    "Drift velocity correction (?)"
    drift_vel::AbstractString

    ".dat/spe file with impurity profile"
    impurity_profile::AbstractString=""

    "offset of detector Z=0 in crystal from seed start"
    offset_in_mm::Real=-1

    "Path to crystal jsons: alternative to .dat file"
    crystal_metadata_path::AbstractString = ""

    "Name for cached simulation. Not caching if empty string"
    cached_name::AbstractString = ""    
end


SiggenSimulator(::Any) = throw(ErrorException("SiggenSimulator requires MJDSigGen to be loaded, e.g. via `import MJDSigGen`"))

# function SSDSimulator(sim_conf::LegendGeSimConfig)
function SSDSimulator(simulation_settings::PropDict)
    if(!haskey(simulation_settings, :crystal_metadata_path))
        simulation_settings[:crystal_metadata_path] = ""
        # simulation_settings.crystal_metadata_path = crystal_metadata_path
        # ToDo: move somewhere else - irrelevant if only geometry is constructed, or when simulation read from cache
        @warn "No crystal metadata path given. Simulation with dummy constant impurity density."
    # elseif !ispath(simulation_settings.crystal_metadata_path)
    #     @error "The path to crystal metadata you provided is not valid! ($(simulation_settings.crystal_metadata_path))"
    else
        @info "Impurity profile information based on $(simulation_settings.crystal_metadata_path)"
    end

    coord = haskey(simulation_settings, :coordinates) ? simulation_settings.coordinates : "cylindrical"
    if !(coord in ["cartesian", "cylindrical"])
        error("$coord coordinates not implemented!\n Available: cartesian, cylindrical")
    end

    comp = haskey(simulation_settings, :computation) ? simulation_settings.computation : "2D"
    if !(comp in ["2D", "3D"])
        error("$comp computation not implemented!\n Available: 2D, 3D")
    end

    time_step = haskey(simulation_settings, :time_step) ? simulation_settings.time_step*u"ns" : 1u"ns"
    diff = haskey(simulation_settings, :diffusion) ? simulation_settings.diffusion : false
    selfrep = haskey(simulation_settings, :self_repulsion) ? simulation_settings.self_repulsion : false
    number_of_carriers = haskey(simulation_settings, :number_of_carriers) ? simulation_settings.number_of_carriers : 1

    SSDSimulator(coord, comp, simulation_settings.crystal_metadata_path, simulation_settings.cached_name,
        time_step, diff, selfrep, number_of_carriers)
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
function simulate_waveforms(stp_events::Table, detector::SolidStateDetectors.Simulation, simulator::SSDSimulator)
    @info("~.~.~ SolidStateDetectors")
    contact_charge_signals = SolidStateDetectors.simulate_waveforms(
            stp_events,
            detector,
            max_nsteps = 20000,
            Δt = simulator.time_step,
            diffusion = simulator.diffusion,
            self_repulsion = simulator.self_repulsion,
            number_of_carriers = simulator.number_of_carriers,
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