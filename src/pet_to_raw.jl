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
    simulation_settings::PropDict, setup_settings::PropDict, noise_settings::PropDict; n_waveforms::Int = 0)
    ## step 1: stepping information (pet->stp)
    # det meta just for geometry
    stp_table = pet_to_stp(detector_metadata, pet_filename)
    stp_of_interest = n_waveforms == 0 ? stp_table : stp_table[1:n_waveforms]

    ## step 2: simulate pulses (stp->pss)
    env = Environment(environment_settings)
    simulator = PSSimulator(simulation_settings)
    noise_model = NoiseModel(noise_settings)

    # pss_table, pss_truth = stp_to_pss(stp_table, det_meta, env, ps_simulator, noise_model, config_name)
    # ToDo: overwrite
    pss_table, pss_truth = stp_to_pss(stp_of_interest, detector_metadata, env, simulator, noise_model)

    ## step 3: simulate DAQ
    elec_chain = ElecChain(setup_settings) # needs preamp and fadc
    trigger = Trigger(setup_settings.trigger)
    daq = DAQ(setup_settings.daq)
    
    pss_to_raw(pss_table, pss_truth, elec_chain, trigger, daq, noise_model)
end


# user launches directly inputting separate dicts
function simulate_raw(detector_metadata::AbstractString, pet_filename::AbstractString, environment_settings::Union{Dict,PropDict},
    simulation_settings::Union{Dict,PropDict}, setup_settings::Union{Dict,PropDict}, noise_settings::Union{Dict,PropDict}=Dict("type"=>"none"); n_waveforms::Int = 0)
    pet_to_raw(readlprops(detector_metadata), pet_filename,
        PropDict(environment_settings), PropDict(simulation_settings), PropDict(setup_settings), PropDict(noise_settings); n_waveforms)    
end

# user launches with all settings in one dict
function simulate_raw(detector_metadata::AbstractString, pet_filename::AbstractString, all_settings::Union{Dict,PropDict}; n_waveforms::Int = 0)
    all_set = PropDict(all_settings)
    noise_settings = haskey(all_set, :noise_model) ? all_set.noise_model : Dict("type" => "none")
    simulate_raw(detector_metadata, pet_filename, all_set.environment, all_set.simulation, all_set.setup, noise_settings; n_waveforms)    
end

# user launches with all settings in json
function simulate_raw(detector_metadata::AbstractString, pet_filename::AbstractString, all_settings::AbstractString; n_waveforms::Int = 0)
    simulate_raw(detector_metadata, pet_filename, readlprops(all_settings); n_waveforms)
end
