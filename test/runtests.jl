# This file is a part of LegendGeSim.jl, licensed under the MIT License (MIT).

using Test

using SolidStateDetectors
using LegendGeSim
using LegendTestData
using Unitful

testdata_path = joinpath(LegendTestData.legend_test_data_path(), "data", "legend", "metadata", "hardware", "detectors", "germanium", "diodes")

test_dict = Dict{String, typeof(1.0u"pF")}(
    "B99000A.json" => -7.13u"pF",
    "V99000A.json" => -3.55u"pF"
)

@testset "Package LegendGeSim" begin
    for (filename, capacitance) in test_dict
        detector_metadata_filename = joinpath(testdata_path, filename)
        sim_settings_ssd_filename = joinpath(dirname(dirname(pathof(LegendGeSim))), "examples/configs/SSD_NoiseSim.json")
        
        @testset "Field Simulation of $filename (SSD)" begin
            config = LegendGeSim.load_config(detector_metadata_filename, sim_settings_ssd_filename);
            sim = LegendGeSim.simulate_fields(config)
            C = LegendGeSim.capacitance_matrix(sim)

            @testset "Capacitances" begin   
                @test isapprox(C[1,2], capacitance, atol = 0.2u"pF")
            end
        end
    end 
end

