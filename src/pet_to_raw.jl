"""
    pet_to_raw(pet_file_path, sim_config, config_name)

AbstractString, PropDict, AbstractString -> Table

[WIP] Full simulation chain pet->stp->pss->raw
    based on simulated energy depositions contained in the HDF5 file
    found in <pet_file_path> and simulation settings given in <sim_config>.

The output name <config_name> is used to construct filenames for cached
    simulation files (currently the same as the simulation config basename)
"""
function pet_to_raw(sim_config::PropDict)
    # function pet_to_raw(pet_file_path::AbstractString, sim_config::PropDict, config_name::AbstractString)
    ## step 1: stepping information
    stp_table = pet_to_stp(sim_config)


    ## step 2: simulate pulses
    # det_meta = PropDicts.read(PropDict, sim_config.detector_metadata)
    # env = Environment(sim_config)
    # ps_simulator = PSSimulator(sim_config)
    # noise_model = NoiseModel(sim_config)

    # pss_table, pss_truth = stp_to_pss(stp_table, det_meta, env, ps_simulator, noise_model, config_name)
    # launch with sim config updated with stp table as new input file
    pss_table, pss_truth = stp_to_pss(load_config(stp_table, sim_config))

    ## step 3: simulate DAQ
    # elec_chain = ElecChain(sim_config)
    # trigger = Trigger(sim_config)
    # daq = DAQ(sim_config)

    # update sim config
    # sim_config = load_config(pss_table, sim_config)
    
    # pss_to_raw(pss_table, pss_truth, elec_chain, trigger, daq, noise_model)
    # launch with sim config updated with pss table as new input file
    pss_to_raw(load_config(pss_table, sim_config), pss_truth)
end



function pet_to_raw(pet_input_fullpath::AbstractString, det_meta_fullpath::AbstractString, sim_config_filename::AbstractString)
    pet_to_raw(load_config(pet_input_fullpath, det_meta_fullpath, sim_config_filename))
end



function pet_to_pss(sim_config::PropDict)
    ## step 1: stepping information
    stp_table = pet_to_stp(sim_config)

    # pss_table, pss_truth = stp_to_pss(stp_table, det_meta, env, ps_simulator, noise_model, config_name)
    # launch with sim config updated with stp table as new input file
    pss_table, pss_truth = stp_to_pss(load_config(stp_table, sim_config))

    pss_table, pss_truth 
end


function pet_to_pss(pet_input_fullpath::AbstractString, det_meta_fullpath::AbstractString, sim_config_filename::AbstractString)
    pet_to_pss(load_config(pet_input_fullpath, det_meta_fullpath, sim_config_filename))
end
