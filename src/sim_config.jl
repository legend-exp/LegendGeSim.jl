function load_config(sim_conf_file::AbstractString)
    @info "Simulation config taken from: $sim_conf_file"
    PropDicts.read(PropDict, sim_conf_file)
end


