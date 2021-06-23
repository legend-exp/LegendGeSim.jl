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
    pss_table, pss_truth = stp_to_pss(sim_config)

    ## step 3: simulate DAQ
    # elec_chain = ElecChain(sim_config)
    # trigger = Trigger(sim_config)
    # daq = DAQ(sim_config)
    
    # pss_to_raw(pss_table, pss_truth, elec_chain, trigger, daq, noise_model)
    pss_to_raw(sim_config)
end


function pet_to_raw(pet_file_path::AbstractString, sim_config_filename::AbstractString)
    # sim_config = load_config(sim_config_filename)    
    sim_config = propdict(sim_config_filename)

    # currently based on simulation config name
    # in principle can be given by user or determined in any other way
    config_name = splitext(basename(sim_config_filename))[1]

    pet_to_raw(pet_file_path, sim_config, config_name)
end



