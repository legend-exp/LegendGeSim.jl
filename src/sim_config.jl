# struct LegendGeSimConfig{DT}
#     dict::DT
# end

# using Base: cache_dependencies
# """
#     load_config(input_file, det_metadata, sim_config_file)

# AbstractString or Table, AbstractString, AbstractString -> LegendGeSimConfig

# Merge simulation inputs and settings into one simulation config to pass around
# """
# function load_config(det_metadata::AbstractString, sim_config_file::AbstractString)
#     # construct default simulation cached name
#     cached_name = filename(sim_config_file)
    
#     # add detector info as simulation input
#     # add default cached name which is the name of the config file itself
#     sim_inputs = PropDict(
#         :detector_metadata => propdict(det_metadata),
#         :simulation => PropDict( :cached_name => cached_name )
#     )

#     # contains all simulation settings (method SSD/siggen, optional cached file name etc.)
#     # except not detector info
#     sim_config = propdict(sim_config_file)    
    
#     # merge
#     # after merging, if cached_name was present in sim_config, it will override the default one above
#     # so, if no cached name is provided SSD will cache anyway -> dangerous cause then modify config but will load old!! ToDo!!
#     sim_total = merge(sim_inputs, sim_config)

#     # workaround because PropDict does not seem to be mutable
#     # add detector name prefix to cached name -> ToDo!! dangerous!
#     temp = PropDict(:simulation => PropDict(:cached_name => "$(filename(det_metadata))_$(sim_total.simulation.cached_name)"))
#     # replace/update cached name
#     sim_total = merge(sim_total, temp)
#     @info "(Future) filename for cached simulation: $(sim_total.simulation.cached_name)"

#     LegendGeSimConfig(sim_total)
# end

# ----------------------------------------------------------------

# ToDo: add units as part of variable type?
@with_kw struct Environment
    "Medium"
    medium::String = ""

    "Crystal temperature in K"
    crystal_temperature::Real = 0

    "Oprating voltage in V"
    operating_voltage::Real = 0
end


"""
Simulation parameters related to the environment
"""
function Environment(env_conf::PropDict)
    Environment(
        env_conf.medium,
        env_conf.crystal_temperature_in_K,
        haskey(env_conf, :operating_voltage_in_V) ? env_conf.operating_voltage_in_V : 0
    )
end