"""
    pet_to_raw(pet_file_path, sim_config, config_name)

AbstractString, PropDict, AbstractString -> Table

[WIP] Full simulation chain pet->stp->pss->raw
    based on simulated energy depositions contained in the HDF5 file
    found in <pet_file_path> and simulation settings given in <sim_config>.

The output name <config_name> is used to construct filenames for cached
    simulation files (currently the same as the simulation config basename)
"""
function pet_to_raw(detector_metadata::PropDict, pet_filename::AbstractString, environment_settings::PropDict, simulation_settings::PropDict, setup_settings::PropDict)
    ## step 1: stepping information (pet->stp)
    # det meta just for geometry
    stp_table = pet_to_stp(detector_metadata, pet_filename)


    ## step 2: simulate pulses (stp->pss)
    env = Environment(PropDict(environment_settings))
    simulator = PSSimulator(PropDict(simulation_settings))

    # det_meta = PropDicts.read(PropDict, sim_config.detector_metadata)
    # env = Environment(sim_config)
    # ps_simulator = PSSimulator(sim_config)

    # pss_table, pss_truth = stp_to_pss(stp_table, det_meta, env, ps_simulator, noise_model, config_name)
    # ToDo: overwrite
    # ToDo: noise model ! ! ! ! !
    pss_table, pss_truth = stp_to_pss(stp_table, detector_metadata, env, simulator)

    ## step 3: simulate DAQ
    setup = PropDict(setup_settings)
    elec_chain = ElecChain(setup) # needs preamp and fadc
    trigger = Trigger(setup.trigger)
    daq = DAQ(setup.daq)
    # noise_model = NoiseModel(setup.noise)
    # TEMP ! !
    noise_model = NoiseFromSim(0)

    # update sim config
    # sim_config = load_config(pss_table, sim_config)
    
    # pss_to_raw(pss_table, pss_truth, elec_chain, trigger, daq, noise_model)
    # launch with sim config updated with pss table as new input file
    # ToDo: does not need to know simulator!!
    pss_to_raw(pss_table, pss_truth, simulator, elec_chain, trigger, daq, noise_model)
end


# provided by user through dicts in notebook
function simulate_raw(detector_metadata::AbstractString, pet_filename::AbstractString, environment_settings::Dict, simulation_settings::Dict, setup_settings::Dict)
    pet_to_raw(propdict(detector_metadata), pet_filename, PropDict(environment_settings), PropDict(simulation_settings), PropDict(setup_settings))    
end

function simulate_raw(detector_metadata::AbstractString, pet_filename::AbstractString, all_settings::Dict)
    all_set = PropDict(all_settings)
    pet_to_raw(detector_metadata, pet_filename, all_set.environment, all_set.simulation, all_set.setup)    
end

# ToDo provided through json


# function pet_to_pss(sim_config::PropDict)
#     ## step 1: stepping information
#     stp_table = pet_to_stp(sim_config)

#     # pss_table, pss_truth = stp_to_pss(stp_table, det_meta, env, ps_simulator, noise_model, config_name)
#     # launch with sim config updated with stp table as new input file
#     pss_table, pss_truth = stp_to_pss(load_config(stp_table, sim_config))

#     pss_table, pss_truth 
# end


# function pet_to_pss(pet_input_fullpath::AbstractString, det_meta_fullpath::AbstractString, sim_config_filename::AbstractString)
#     pet_to_pss(load_config(pet_input_fullpath, det_meta_fullpath, sim_config_filename))
# end
