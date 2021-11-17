# This file is a part of LegendGeSim.jl, licensed under the MIT License (MIT).

using Test

using SolidStateDetectors
using LegendGeSim
using LegendTestData
using Unitful

testdata_path = joinpath(LegendTestData.legend_test_data_path(), "data", "ldsim")

@testset "Package LegendGeSim" begin
    detector_metadata_file = joinpath(testdata_path, "invcoax-metadata.json") # path to my detector metadata
    simulation_config_file = joinpath(dirname(dirname(pathof(LegendGeSim))), "examples/configs/SSD_NoiseSim.json")
    fields_sim_file = "det_fields_ssd" # file prefix for automated storing of simulated fields

    simulation_settings = LegendGeSim.SSDSimulator(LegendGeSim.propdict(simulation_config_file))
    environment_settings = LegendGeSim.Environment(LegendGeSim.propdict(simulation_config_file))
    
    @testset "Field Simulation (SSD)" begin
        sim = LegendGeSim.simulate_fields(detector_metadata_file, environment_settings, simulation_settings, fields_sim_file; overwrite = true)

        C = calculate_capacitance_matrix(sim) 
        active_volume = get_active_volume(sim.point_types)

        @testset "Capacitances" begin   
            @test isapprox(ustrip(active_volume), 230.5, atol = 1.0)
            @test isapprox(ustrip(C[1,2]), -3.25, atol = 0.2)
        end
    end
end 

