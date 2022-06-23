"""
The DAQ supertype corresponds to the component of the 
    DAQ and electronics setup that saves the waveform after 
    it triggered
"""
abstract type DAQ end

"""
GenericDAQ is a dummy DAQ model that stores the waveform
    after it triggered, saving <baseline_length> samples before the
    trigger index, and in total <nsamples> samples from start to end
"""
@with_kw struct GenericDAQ <: DAQ
    "Number of samples i.e. length of the final stored waveform"
    nsamples::Int = 4000

    "Length of stored baseline in number of samples"
    baseline_length::Int = 400
end


"""
    GenericDAQ(sim_conf)

PropDict -> GenericDAQ 

Construct a GenericDAQ struct based on given simulation configuration 
"""
function GenericDAQ(sim_conf::PropDict)
    GenericDAQ(
        nsamples = sim_conf.setup.daq.nsamples,
        baseline_length = sim_conf.setup.daq.baseline_length
    )
end


"""
    DAQ(sim_conf)

PropDict -> <DAQ>

Construct a DAQ supertype struct based on given simulation configuration.
Type of returned instance depends on settings in <sim_conf>
Currently only one type of DAQ available (GenericDAQ),
    rendering this function temporarily redundant.
"""
function DAQ(sim_config::PropDict)
    GenericDAQ(sim_config)
end

# -------------------------------------------------------------------

"""
    simulate(wf, trigger_index, daq)

RDWaveform, Int, DAQ -> RDWaveform

Simulate how the DAQ stores the waveform after it receives a trigger.
"""
function simulate(wf::RDWaveform, trigger_index::Int, daq::GenericDAQ)
    T = Float32 # This should be somehow defined and be passed properly
    # the stored waveform will start from zero
    # and go until the number of samples DAQ is configured to save
    ts = range(T(0)u"ns", step = step(wf.time), length = daq.nsamples)
    
    # in case the waveform didn't trigger (which is what 0 index stands for)...
    if(trigger_index == 0)
        # ...return all zeros
        return RDWaveform(ts, similar(wf.signal[1:daq.nsamples])) # zeros
    else
        # ...otherwise store the waveform
        iStart = trigger_index - daq.baseline_length
        return RDWaveform(ts, wf.signal[iStart : iStart + daq.nsamples-1])
    end
end
