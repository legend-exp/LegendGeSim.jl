"""
    pss_to_raw(pss_table, pss_truth, simulation_settings, elec_chain, trigger, daq, noise_model)

Table, Table, SSDSimulator, ElecChain, Trigger, DAQ, NoiseModel -> Table

Simulate effects of the electronics chain <elec_chain>, <trigger> and <daq> settings
    on the waveforms contained in <pss_table>.
Noise is simulated based on the given <noise_model>.
Construct a table with the format identical to data raw tier.
Currently timing information in <pss_truth> is used for dummy timestamps in the output
    raw tier table.
"""
function pss_to_raw(pss_table::Table, pss_truth::Table, elec_chain::ElecChain, trigger::Trigger, daq::DAQ, noise_model::NoiseModel)

    @info "---------------------- pss -> raw (realistic waveforms)"

    result = process_waveforms(pss_table, elec_chain, trigger, daq, noise_model)
   
    ## construct array of waveforms for hdf5
    wf_final = ArrayOfRDWaveforms(result.wf_array)
    # why am I doing this?
    wf_final = ArrayOfRDWaveforms((wf_final.time, VectorOfSimilarVectors(wf_final.signal)))

    raw_table = Table(
        baseline = result.baseline,
        channel = pss_table.channel,
        energy = result.online_energy,
        ievt = pss_table.ievt,
        # why am i using baseline here? could be any column?
        numtraces = ones(length(result.baseline)), # number of triggered detectors (1 for HADES)
        packet_id = zeros(length(result.baseline)), # means to packet losses
        timestamp = getindex.(pss_truth.thit, 1), # frist MC truth hit time of each event?
        tracelist = VectorOfVectors([[1] for idx in 1:length(result.baseline)]), # lists of ADCs that triggered, 1 for HADES all the time
        waveform = wf_final,
        wf_max = maximum.(wf_final.signal),
        # I was told that std is std of the whole waveform, not just the baseline
        wf_std = std.(wf_final.signal)
    )

    # probably not needed because one can look this up in the simulation config file
    # TODO
    # electruth = Table(
    #     # fill in noise and DAQ parameters
    # )

    raw_table
end

# wrapper for user reading from pre-saved pss file
function pss_to_raw(pss_file::AbstractString, setup_settings::Dict, noise_settings::Dict=Dict("type"=>"none"))
    pss_h5 = h5open(pss_file, "r")
    pss_table = LegendHDF5IO.readdata(pss_h5, "pss/pss")
    pss_truth = LegendHDF5IO.readdata(pss_h5, "pss/truth")
    close(pss_h5)

    setup = PropDict(setup_settings)
    elec_chain = ElecChain(setup) # needs preamp and fadc
    trigger = Trigger(setup.trigger)
    daq = DAQ(setup.daq)
    noise_model = NoiseModel(PropDict(noise_settings))
    # TEMP ! !
    # noise_model = NoiseFromSim(0u"keV")

    pss_to_raw(pss_table, pss_truth, elec_chain, trigger, daq, noise_model)
end


# ------------------------------





"""
    process_waveforms(pss_table, elec_chain, trigger, daq, noise_model) 

Table, GenericElecChain, Trigger, DAQ, NoiseModel -> Table  

Simulate effects of the electronics chain <elec_chain>, <trigger> and <daq> settings
    on the waveforms contained in <pss_table>.
Noise is simulated based on the given <noise_model>.
The output table contains the resulting DAQ waveforms, simulated online energy after tigger,
    baselines and their RMS, needed for the raw tier table.
"""
function process_waveforms(pss_table::Table, elec_chain::GenericElecChain, trigger::Trigger, daq::DAQ, noise_model::NoiseModel)
    T = Float32 # This should be somehow defined and be passed properly
    # Create arrays to be filled with results, and online energy
    n_waveforms = size(pss_table.waveform, 1)
    result = Table(
        wf_array = Array{RDWaveform}(undef, n_waveforms),
        online_energy = Array{typeof(1.0 * energy_unit)}(undef, n_waveforms),
        baseline = Array{T}(undef, n_waveforms),
        baseline_rms = Array{T}(undef, n_waveforms)
    )

    ## ---- temporary

    # in the future won't be needed, because will require trigger threshold in input ?
    # right now to make things easier, I calculate it based on noise levels

    # need to think about what to do with preamp gain if baselines are added to simulated wf 
    # right now calculating gain based on offset in data baseline 
    # has to be before trigger threshold obv!
    elec_chain.preamp.gain = preamp_gain(elec_chain.preamp, noise_model)
    trigger.threshold = trigger_threshold(trigger, elec_chain.preamp, noise_model)

    @info "Processing waveforms..."
    for i = 1:n_waveforms
        # for some reason ProgressMeter doesn't work here
        if (i % 100 == 0)
            println("$i / $n_waveforms")
        end

        wf = pss_table.waveform[i]

        # all the operations below assume unitless wf
        # currently stp->pss for siggen returns unitless (?), SSD in eV
        # strip, ToDo: figure out what to do here
        wf = RDWaveform(wf.time, ustrip.(wf.signal))

        ## invert the pulse if it came from the n+ contact
        sign = wf.signal[end] < 0 ? -1 : 1
        wf = RDWaveform(wf.time, sign * wf.signal)

        ## negative values control - TEMP
        # for now simply replace with zero 
        wf = RDWaveform(wf.time, remove_negative.(wf.signal))

        # extend tail and baseline long enough for future DAQ processing
        wf = add_tail_and_baseline(wf, elec_chain.fadc, daq)

        # obtain current from charge
        wf = differentiate(wf)

        # simulate electronics chain
        wf = simulate(wf, elec_chain, noise_model)

        # simulate trigger
        trigger_index, online_energy = simulate(wf, trigger)

        # simulate DAQ
        wf = simulate(wf, trigger_index, daq)

        ## fill the wf array
        result.wf_array[i] = wf

        # if online energy is zero, means didn't trigger
        # convert online energy to keV
        result.online_energy[i] = trigger_index == 0 ? 0u"keV" : uconvert(u"keV", online_energy / elec_chain.preamp.gain * germanium_ionization_energy)
        result.baseline[i], result.baseline_rms[i] = mean_and_std(result.wf_array[i].signal[1:daq.baseline_length])
    end

    result
end