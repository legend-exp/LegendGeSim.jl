function load_config(sim_conf_file::AbstractString)
    @info "Simulation config taken from: $sim_conf_file"
    PropDicts.read(PropDict, sim_conf_file)
end


@with_kw struct Environment
    crystal_t::Real = 0
    op_voltage::Real = 0
end


"""
Simulation parameters related to the environment
"""
function Environment(sim_conf::PropDict)
    Environment(sim_conf.environment.crystal_t, sim_conf.environment.op_voltage)
end