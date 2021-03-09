function load_config(sim_conf_file::AbstractString)
    @info "Simulation config taken from: $sim_conf_file"
    sim_conf = PropDicts.read(PropDict, sim_conf_file)

    sim_conf
end

# sim_config = load_config()

