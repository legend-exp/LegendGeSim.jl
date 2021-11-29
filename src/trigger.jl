"""
The Trigger supertype corresponds to the trigger component in DAQ.
"""
abstract type Trigger end


"""
The TrapFilter trigger is a dummy simulation of the trapezoidal filter
    algorithm to produce a trigger if the waveform amplitude passes
    a certain threshold.

This struct is currently mutable because in case of NoiseFromData modelling,
    we need to calculate the threshold post-factum after analysing the noise
    levels in data.    
"""
@with_kw mutable struct TrapFilter <: Trigger
    "Lenghs of the three trap filter windows in number of samples"
    window_lengths::NTuple{3, Int} = (250,250,250)

    "Trigger threshold in keV"
    threshold_keV::typeof(1.0*energy_unit) = 0u"keV"

    "Trigger threshold after conversion"
    threshold::Real = 0

end


"""
    TrapFilter(sim_conf)

PropDict -> TrapFilter 

Construct a TrapFilter instance based on simulation configuration given in <sim_conf>.
"""
function TrapFilter(sim_conf::PropDict)
    T = Float32 # This should be somehow defined and be passed properly
    TrapFilter(
        window_lengths = (sim_conf.setup.trigger.window_lengths[1],
            sim_conf.setup.trigger.window_lengths[2], sim_conf.setup.trigger.window_lengths[3]),
        threshold_keV = T(sim_conf.setup.trigger.threshold)u"keV"
    )
end


"""
    Trigger(sim_conf)

PropDict -> <Trigger>    

Construct a Trigger supertype instance based on settings given in <sim_conf>.
The returned type depends on the given settings.

Currently only TrapFilter type is implemented.

"""
function Trigger(sim_conf::PropDict)
    if sim_conf.setup.trigger.type == "trapezoidal"
        TrapFilter(sim_conf)
    else
        @info "Trigger type $(sim_config.setup.trigger.type) not implemented!\n
        Available type: trapezoidal"
    end
end

# -------------------------------------------------------------------

"""
    trap_filter(wf_values, sample_idx, trap_filter)

AbstractVector, Int, TrapFilter -> Real, Bool 

Simulate the output of the trapezoidal filter algorithm
    applied to <wf_values> when the first sliding window starts
    at <sample_idx>.
The trapezoidal algorithm parameters are contained in <trap_filter.
The returned values are the online energy estimate corresponding to the
    difference between the 1st and the 3rd sliding window, and a Bool
    correspondong to whether there is a trigger issued or not.    
"""
function trap_filter(wf_values::AbstractVector, sample_idx::Int, trap_filter::TrapFilter)
    # first window
    r1 = sample_idx : sample_idx + trap_filter.window_lengths[1]
    # third window
    r2 = sample_idx + trap_filter.window_lengths[1] + trap_filter.window_lengths[2] : sample_idx + sum(trap_filter.window_lengths)
    # difference between the third and first window
    r = mean(wf_values[r2]) - mean(wf_values[r1])
    r, r >= trap_filter.threshold
end


function trap_filter(wf::RDWaveform, sample_idx::Int, trigger::TrapFilter)
    trap_filter(wf.value, sample_idx, trigger)
end


"""
    simulate(wf, trap_filter)

RDWaveform, TrapFilter -> Int, Real 

Simulate the results of applying <trap_filter> to <wf>
Returns the index on which <wf> triggered, and the estimated
    online energy based on <trap_filter> output.
"""
function simulate(wf::RDWaveform, trigger::TrapFilter)
    T = Float32 # This should be somehow defined and be passed properly
    trigger_window_length = sum(trigger.window_lengths)

    online_filter_output = zeros(T, length(wf.value) - trigger_window_length)
    t0_idx::Int = 0
    trig = false
    
    # replace by a while loop?
    for i in eachindex(online_filter_output)
        online_filter_output[i], trig = trap_filter(wf, i, trigger)
        if trig && t0_idx == 0
            t0_idx = i + trigger.window_lengths[1] + trigger.window_lengths[2]
        end
    end

    t0_idx, maximum(online_filter_output)
end


