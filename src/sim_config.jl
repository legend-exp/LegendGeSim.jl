# using Base: cache_dependencies
"""
    load_config(input_file, det_metadata, sim_config_file)

AbstractString or Table, AbstractString, AbstractString -> PropDict

Merge simulation inputs and settings into one simulation config to pass around
"""
function load_config(det_metadata::AbstractString, sim_config_file::AbstractString)
    # function load_config(input_file::Union{AbstractString, Table}, det_metadata::AbstractString, sim_config_file::AbstractString)
    # construct default cached name
    cached_name = filename(sim_config_file)

    # add detector info as simulation input
    sim_inputs = PropDict(
        :detector_metadata => det_metadata,
        :simulation => PropDict( :cached_name => cached_name )
    )

    # contains only simulation settings
    sim_config = propdict(sim_config_file)    
    
    # merge
    # after merging, if cached_name was present in sim_config, it will override the default one above
    sim_total = merge(sim_inputs, sim_config)

    # add detector name to cached filename
    # what? this doesn't work?? can i make it mutable?
    # sim_total.pss.cached_name = "$(filename(sim_total.detector_metadata))_$(sim_total.pss.cached_name)"
    # workaround
    temp = PropDict(:pss => PropDict(:cached_name => "$(filename(sim_total.detector_metadata))_$(sim_total.simulation.cached_name)"))
    sim_total = merge(sim_total, temp)
    @info "(Future) cached basename: $(sim_total.simulation.cached_name)"

    # load_config(input_file, sim_total)
    sim_total
end


function load_config(input_file::Union{AbstractString, Table}, sim_config::PropDict)
    # add input
    sim_inputs = PropDict(
        :input_file => input_file, 
    )

    # merge
    # new input file will override previous input if present
    sim_total = merge(sim_config, sim_inputs)

    #@info "Your final simulation configuration" sim_total
    
    sim_total
end


function load_config(input_file::Union{AbstractString, Table}, det_metadata::AbstractString, sim_config_file::AbstractString)
    sim_config = load_config(det_metadata, sim_config_file)

    load_config(input_file, sim_config)
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
function Environment(sim_conf::PropDict)
    Environment(
        sim_conf.environment.medium,
        sim_conf.environment.crystal_t,
        sim_conf.environment.op_voltage
    )
end
