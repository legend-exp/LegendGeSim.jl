struct LegendGeSimConfig{DT}
    dict::DT
end

# using Base: cache_dependencies
"""
    load_config(input_file, det_metadata, sim_config_file)

AbstractString or Table, AbstractString, AbstractString -> LegendGeSimConfig

Merge simulation inputs and settings into one simulation config to pass around
"""
function load_config(det_metadata::AbstractString, sim_config_file::AbstractString)
    # construct default simulation cached name
    cached_name = filename(sim_config_file)
    
    # add detector info as simulation input
    sim_inputs = PropDict(
        :detector_metadata => propdict(det_metadata),
        :simulation => PropDict( :cached_name => cached_name )
    )

    # contains only simulation settings
    sim_config = propdict(sim_config_file)    
    
    # merge
    # after merging, if cached_name was present in sim_config, it will override the default one above
    sim_total = merge(sim_inputs, sim_config)

    # workaround because PropDict does not seem to be mutable
    # add detector name prefix to cached name
    temp = PropDict(:simulation => PropDict(:cached_name => "$(filename(det_metadata))_$(sim_total.simulation.cached_name)"))
    # replace/update cached name
    sim_total = merge(sim_total, temp)
    @info "(Future) filename for cached simulation: $(sim_total.simulation.cached_name)"

    LegendGeSimConfig(sim_total)
end


@with_kw struct Environment
    "Medium"
    medium::String = "vacuum"

    "Crystal temperature in K"
    crystal_t::Real = 0

    "Oprating voltage in V"
    op_voltage::Real = 0
end


"""
Simulation parameters related to the environment
"""
function Environment(sim_conf::LegendGeSimConfig)
    Environment(
        sim_conf.dict.environment.medium,
        sim_conf.dict.environment.crystal_t,
        sim_conf.dict.environment.op_voltage
    )
end
