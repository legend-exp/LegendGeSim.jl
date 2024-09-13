# ToDo: add units as part of variable type?
@with_kw struct Environment
    "Medium"
    medium::String = ""

    "Crystal temperature in K"
    crystal_temperature::Real = 0

    "Oprating voltage in V"
    operating_voltage::Real = 0
    
    "Quickfix DL"
    dl::Union{Real,AbstractString} = 0
end


"""
Simulation parameters related to the environment
"""
function Environment(env_conf::PropDict)
    Environment(
        env_conf.medium,
        env_conf.crystal_temperature_in_K,
        haskey(env_conf, :operating_voltage_in_V) ? env_conf.operating_voltage_in_V : 0,
        haskey(env_conf, :dl) ? env_conf.dl : 0
    )
end
