# Random.seed!(123) # only for testing


"""
    mstp_to_mcpss(mc_events, sim_config_file)

Simulate waveforms based on given MC events for a given configuration.
Returns the resulting simulated waveforms and MC truth

mc_events: table in the mcstp format (output of g4_to_mcstp)
sim_config_file: simulation configuration json file
"""
function mcstp_to_mcpss(mc_events::Table, sim_config_file::AbstractString)
    sim_config = LegendGeSim.load_config(sim_config_file)
    mcstp_to_mcpss(mc_events, sim_config)

end

function mcstp_to_mcpss(mc_events::Table, sim_config::PropDict)
    ps_simulator = PSSimulator(sim_config)
    @info("Reading geometry for $(sim_config.detector)")
    det_config = detector_config(sim_config.detector, ps_simulator)
    noise_model = NoiseModel(sim_config)

    mcstp_to_mcpss(mc_events, det_config, ps_simulator, noise_model)

end


"""
    mstp_to_mcpss(mc_events, sim_config)

Simulate waveforms based on given MC events for a given configuration.
Returns the resulting simulated waveforms and MC truth

mc_events: table in the mcstp format (output of g4_to_mcstp)
sim_config: PropDict object with simulation configuration
"""
function mcstp_to_mcpss(mc_events::Table, detector_config::Dict, ps_simulator::PSSimulator, noise_model::NoiseModel)
    # add fano noise, don't add if data noise is applied later
    # should we do this before siggen as well? or does siggen do it by itself?
    # noise_model = NoiseModel(sim_config)
    mc_events = fano_noise(mc_events, detector_config, noise_model)

    mcpss_table, mcpss_mctruth = simulate_wf(mc_events, detector_config, ps_simulator)

    mcpss_table, mcpss_mctruth
end


# function detector_simulation(det_path::AbstractString, det_name::AbstractString, ::SiggenSimulator)
#     PropDicts.read(PropDict, joinpath(det_path, det_name*".json"))
# end


# """
#     detector_simulation(det_path, det_name)

# Read cached h5 detector simulation if exists,
# otherwise simulate detector based on json geometry.

# det_path: path to detector json file
# det_name: detector name without extension (e.g. "IC160A")

# Output: SSD detector simulation object
# """
# function detector_simulation(detector::AbstractString, ::SSDSimulator)

#     det_h5 = joinpath(det_path, det_name*".h5f")
#     if isfile(det_h5)
#         @info "Reading $det_name simulation from cached h5"
#         simulation = SolidStateDetectors.ssd_read(det_h5, Simulation)
#     else
#         @info "Simulating $det_name from scratch"
#         simulation = simulate_detector(det_path, det_name)
#     end

#     simulation
# end





