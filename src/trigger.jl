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

LegendGeSimConfig -> TrapFilter 

Construct a TrapFilter instance based on simulation configuration given in <sim_conf>.
"""
function TrapFilter(trigger_settings::PropDict)
    T = Float32 # This should be somehow defined and be passed properly
    TrapFilter(
        window_lengths = (trigger_settings.window_lengths[1],
            trigger_settings.window_lengths[2],
            trigger_settings.window_lengths[3]),
        threshold_keV = T(trigger_settings.threshold)u"keV"
    )
end


"""
    Trigger(sim_conf)

LegendGeSimConfig -> <Trigger>    

Construct a Trigger supertype instance based on settings given in <sim_conf>.
The returned type depends on the given settings.

Currently only TrapFilter type is implemented.

"""
function Trigger(trigger_settings::PropDict)
    if trigger_settings.type == "trapezoidal"
        TrapFilter(trigger_settings)
    else
        error("Trigger type $(trigger_settings.type) not implemented!\n
        Available type: trapezoidal")
    end
end

# -------------------------------------------------------------------

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
    
    fi = fltinstance(TrapezoidalChargeFilter{T}(reverse(trigger.window_lengths)...), smplinfo(wf)) 
    online_filter_output = similar(wf.signal, length(wf.signal) - trigger_window_length + 1)
    rdfilt!(online_filter_output, fi, wf.signal)
    
    t0idx = findfirst(online_filter_output .>= trigger.threshold) + fi.ngap + fi.navg2 - 1
    t0idx, T(maximum(online_filter_output))
end


