
"""
    add_tail_and_baseline(wf, fadc, daq)

RDWaveform, GenericFADC, GenericDAQ -> RDWaveform  

Extend tail and add baseline of <wf> based on <fadc> and <daq> parameters.

Extended conservatively to ensure enough samples for the future
    sliding trap filter window and the resulting DAQ baseline.
"""
function add_tail_and_baseline(wf::RDWaveform, fadc::GenericFADC, daq::GenericDAQ)
    factor::Integer = uconvert(u"ns", fadc.Δt) / uconvert(u"ns", step(wf.time))
    # double the number of samples should be enough
    # ideally infinite flat baseline and tail
    SolidStateDetectors.add_baseline_and_extend_tail(wf, daq.baseline_length*factor*2, daq.nsamples*factor*2)
end


"""
    differentiate(wf)

RDWaveform -> RDWaveform    

Differentiate a waveform using a Biquad filter
"""
function differentiate(wf::RDWaveform)
    diff_biquad_filter = RadiationDetectorDSP.differentiator_filter(T(1)) # We want a gain of 1 here
    filter_output = filt(diff_biquad_filter, wf.value)
    # Implement a method do parse a RDWaveform to `filt`? (ToDo)
    RDWaveform(wf.time, filter_output)
end