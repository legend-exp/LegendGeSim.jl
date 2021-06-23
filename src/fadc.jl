"""
The FADC supertype corresponds to the FADC component
    of the electronics chain in the DAQ setup.
"""
abstract type FADC end

"""
GenericFADC is a dummy FADC model that performs sampling
    of the input waveform based on the sampling interval <Δt>
"""
@with_kw struct GenericFADC <: FADC
    "Sampling interval of FADC, inverse of sampling rate"
    Δt::typeof(1.0*ns_unit) = uconvert(u"ns", inv(250u"MHz"))

    # to add: number of bits
end


"""
FADC with Flashcam algorithm (ToDo)
"""
struct Flashcam <: FADC end

"""
FADC with Struck algorithm (ToDo)
"""
struct Struck <: FADC end


"""
    GenericFADC(sim_conf)

PropDict -> GenericFADC

Construct a GenericFADC instance based on simulation
    configuration <sim_conf>
"""
function GenericFADC(sim_conf::PropDict)
    GenericFADC(
        Δt = haskey(sim_conf.setup.fadc, :sampling_rate) ? uconvert(u"ns", inv(T(sim_conf.setup.fadc.sampling_rate)u"MHz")) : T(sim_conf.setup.fadc.sampling_interval)u"ns",
    )
end


"""
    FADC(sim_conf)

PropDict -> <FADC>

Construct an FADC supertype instance based on simulation
    configuration <sim_conf>.
Type of <FADC> depends on the type specified in <sim_conf>
    (e.g. generic, Flashcam, Struck)
"""
function FADC(sim_conf::PropDict)
    if sim_conf.setup.fadc.type == "generic"
        GenericFADC(sim_conf)
    else
        @info "FADC type $(sim_conf.setup.fadc.type) not implemented!\n
        Available types: generic\n
        Planned in the future: Flashcam, Struck"              
    end
end

# -------------------------------------------------------------------

"""
    simulate(wf, fadc)   

RDWaveform, GenericFADC -> RDWaveform

Simulate effects of FADC module <fadc> on the waveform <wf>
    (sampling and digitization)
"""
function simulate(wf::RDWaveform, fadc::GenericFADC)   
    # resample -> currently simply pick every Nth sample, more elaborate simulation later
    wf_sampled = wf.value[begin : Int(fadc.Δt/step(wf.time)) : end]
    t_sampled = range(T(0)u"ns", step = fadc.Δt, length = length(wf_sampled))
    wf_daq = RDWaveform(t_sampled, wf_sampled)

    # digitize
    wf_daq = RDWaveform(wf_daq.time, UInt16.(round.(wf_daq.value, digits = 0)))

    wf_daq
end