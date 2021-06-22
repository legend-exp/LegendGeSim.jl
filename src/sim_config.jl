"""
    propdict(json_file)

AbstractString -> PropDict 

Construct a PropDict based on given <json_file>.

1) Why can't I name it function PropDict()?
2) Why doesn't this already exist in PropDicts?
I find the PropDicts.read(PropDict, String) format kinda cumbersome
"""
function propdict(json_file::AbstractString)
    PropDicts.read(PropDict, json_file)
end


@with_kw struct Environment
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
        sim_conf.environment.crystal_t,
        sim_conf.environment.op_voltage
    )
end