"""
    stp_to_pss(stp_table, det_meta, env, ps_simulator, config_name)

Table, PropDict, Env, PSSimulator, NoiseModel, AbstractString -> Table, Table    

* NOTE * Lukas removed NoiseModel from here for some reason (search for fano_noise in old files)

Simulate waveforms based on stepping info given in <stp_table>, 
    LEGEND detector metadata <det_meta>, environment settings given in <env>,
    pulse shape simulation method <ps_simulator>, and <noise_model>

The output is a table with simulated pulses, and a table with simulation truth
    (may be abolished in the future since basically corresponds to stp table)
"""
function stp_to_pss(stp_table::Table, det_meta::PropDict, env::Environment, simulator::PSSimulator, noise_model::NoiseModel)
    @info "---------------------- stp -> pss (ideal pulses)"

    # @info "//\\//\\//\\ Fano noise"
    stp_table = fano_noise(stp_table, det_meta, env, noise_model)

    @info "_||_||_||_ Simulate detector" 
    sim = simulate_detector(det_meta, env, simulator; overwrite = false)

    @info "~.~.~.~.~ Simulate charge pulses"
    simulate_waveforms(stp_table, sim, simulator)
end


# function stp_to_pss(stp_filepath::AbstractString, det_meta_fullpath::AbstractString, sim_config_filename::AbstractString)
#     stp_table = h5open(stp_filepath, "r") do input
#         LegendHDF5IO.readdata(input, "stp")
#     end
    
#     stp_to_pss(stp_table, det_meta_fullpath, sim_config_filename)
# end

# function stp_to_pss(stp_table::Table, sim::Union{Simulation, SigGenSetup}, sim_settings::PSSimulator)
#     # Right now, only use `sim_settings` to dispatch for SSD 
#     # Later, waveform simulation settings could be set in the configuration file and be parsed to `sim_settings`
#     @info "~.~.~.~.~ Simulate charge pulses"
#     simulate_waveforms(stp_table, sim)
# end

# pet->pss
# user launches directly inputting separate dicts
function simulate_pulses(detector_metadata::AbstractString, pet_filename::AbstractString, environment_settings::Union{Dict,PropDict},
    simulation_settings::Union{Dict,PropDict}, noise_settings::Union{Dict,PropDict} = Dict("type" => "none"); n_waveforms::Int = 0)
    ## step 1: stepping information
    det_meta = propdict(detector_metadata)

    stp_table = pet_to_stp(det_meta, pet_filename)
    # truncate if asked by user
    stp_of_interest = n_waveforms == 0 ? stp_table : stp_table[1:n_waveforms]
    # woudl this work? stp_table[1 : n_waveforms == 0 ? end : n_waveforms]

    env = Environment(PropDict(environment_settings))
    simulator = PSSimulator(PropDict(simulation_settings))
    noise_model = NoiseModel(PropDict(noise_settings))

    pss_table, pss_truth = stp_to_pss(stp_of_interest, det_meta, env, simulator, noise_model)

    pss_table, pss_truth 
end

# user launches with all settings in one dict
function simulate_pulses(detector_metadata::AbstractString, pet_filename::AbstractString, all_settings::Union{Dict,PropDict}; n_waveforms::Int = 0)
    all_settings = PropDict(all_settings)
    noise_settings = haskey(all_settings, :noise_model) ? all_settings.noise_model : Dict("type" => "none")
    simulate_pulses(detector_metadata, pet_filename, all_settings.environment, all_settings.simulation, noise_settings; n_waveforms)
end

# user launches with all settings in json
function simulate_pulses(detector_metadata::AbstractString, pet_filename::AbstractString, all_settings::AbstractString; n_waveforms::Int = 0)
    simulate_pulses(detector_metadata, pet_filename, propdict(all_settings); n_waveforms)
end


