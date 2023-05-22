"""
    pet_to_raw(pet_file_path, sim_config, config_name)

AbstractString, PropDict, AbstractString -> Table

[WIP] Full simulation chain pet->stp->pss->raw
    based on simulated energy depositions contained in the HDF5 file
    found in <pet_file_path> and simulation settings given in <sim_config>.

The output name <config_name> is used to construct filenames for cached
    simulation files (currently the same as the simulation config basename)
"""
function pet_to_raw(detector_metadata::PropDict, pet_filename::AbstractString, environment_settings::PropDict,
    simulation_settings::PropDict, setup_settings::PropDict; n_waveforms::Int = 0)
    ## step 1: stepping information (pet->stp)
    # det meta just for geometry
    stp_table = pet_to_stp(detector_metadata, pet_filename)


    ## step 2: simulate pulses (stp->pss)
    env = Environment(PropDict(environment_settings))
    simulator = PSSimulator(PropDict(simulation_settings))

    # pss_table, pss_truth = stp_to_pss(stp_table, det_meta, env, ps_simulator, noise_model, config_name)
    # ToDo: overwrite
    # ToDo: noise model ! ! ! ! !
    pss_table, pss_truth = stp_to_pss(stp_table, detector_metadata, env, simulator; n_waveforms)

    ## step 3: simulate DAQ
    setup = PropDict(setup_settings)
    elec_chain = ElecChain(setup) # needs preamp and fadc
    trigger = Trigger(setup.trigger)
    daq = DAQ(setup.daq)
    # noise_model = NoiseModel(setup.noise)
    # TEMP ! !
    noise_model = NoiseFromSim(0)
    
    # pss_to_raw(pss_table, pss_truth, elec_chain, trigger, daq, noise_model)
    # launch with sim config updated with pss table as new input file
    # ToDo: does not need to know simulator!!
    pss_to_raw(pss_table, pss_truth, elec_chain, trigger, daq, noise_model)
end


# user launches directly inputting separate dicts
function simulate_raw(detector_metadata::AbstractString, pet_filename::AbstractString, environment_settings::Union{Dict,PropDict},
    simulation_settings::Union{Dict,PropDict}, setup_settings::Union{Dict,PropDict}; n_waveforms::Int = 0)
    pet_to_raw(propdict(detector_metadata), pet_filename,
        PropDict(environment_settings), PropDict(simulation_settings), PropDict(setup_settings); n_waveforms=n_waveforms)    
end

# user launches with all settings in one dict
function simulate_raw(detector_metadata::AbstractString, pet_filename::AbstractString, all_settings::Union{Dict,PropDict}; n_waveforms::Int = 0)
    all_set = PropDict(all_settings)
    simulate_raw(detector_metadata, pet_filename, all_set.environment, all_set.simulation, all_set.setup; n_waveforms)    
end

# user launches with all settings in json
function simulate_raw(detector_metadata::AbstractString, pet_filename::AbstractString, all_settings::AbstractString; n_waveforms::Int = 0)
    simulate_raw(detector_metadata, pet_filename, propdict(all_settings); n_waveforms=n_waveforms)
end
