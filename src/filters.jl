"""
    dspjl_differentiator_filter(gain::Real)

Differentiator filter (later to be part of LegendDSP)

The waveform is in units of induced charge (integrated form).
But we want to transform it into its differential form (current).
We will use a BiQuad filter for this: https://en.wikipedia.org/wiki/Digital_biquad_filter

gain: ?
Output: ?
"""
function dspjl_differentiator_filter(gain::Real)
    T = float(typeof(gain))
    Biquad(T(gain), T(-gain), T(0), T(0), T(0))
end



"""
    dspjl_simple_csa_response_filter(τ_rise, τ_decay, gain)

Simulate CSA response using the RC and integrator filters
τ_rise: ?
τ_decay: ?
gain: ?

Output: ?
"""
function dspjl_simple_csa_response_filter(τ_rise::Real, τ_decay::Real, gain::Real = one(τ_rise))
    # TODO: Use a single biquad filter
    T = float(promote_type(promote_type(typeof(τ_rise), typeof(τ_decay)), typeof(gain)))
    dspjl_rc_filter(T(τ_rise)) * dspjl_integrator_cr_filter(T(gain), T(τ_decay))
end


"""
    dspjl_rc_filter(RC)

An `RC` filter made with BiQuad filter

Differentiate a waveform using Biquad filter
(see function dspjl_differentiator_filter)

RC: ?
Output: ?
"""
function dspjl_rc_filter(RC::Real)
    T = float(typeof(RC))
    α = 1 / (1 + RC)
    Biquad(T(α), T(0), T(0), T(α - 1), T(0))
end


"""
    dspjl_integrator_cr_filter(gain, RC)

An `integrator` filter (the inverser of the `dspjl_differentiator_filter`)

gain: ?
RC: ?

Output: ?
"""
function dspjl_integrator_cr_filter(gain::Real, RC::Real)
    T = float(promote_type(typeof(gain), typeof(RC)))
    α = 1 / (1 + RC)
    Biquad(T(gain), T(-α), T(0), T(α - 1), T(0))
end
