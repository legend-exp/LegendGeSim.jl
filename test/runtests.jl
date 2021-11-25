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

    
    @testset "Field Simulation (SSD)" begin
        config = LegendGeSim.load_config(detector_metadata_filename, sim_settings_ssd_filename);

        sim = LegendGeSim.simulate_fields(config)

        C = calculate_capacitance_matrix(sim) 

        @testset "Capacitances" begin   
            @test isapprox(ustrip(C[1,2]), -3.25, atol = 0.2)
        end
    end
end 

