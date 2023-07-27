# This file is a part of LegendGeSim.jl, licensed under the MIT License (MIT).

using Test

using LegendGeSim
using LegendHDF5IO, HDF5
using LegendTestData
using Unitful

# Temporary quick fix to compare RDWaveforms
using LegendGeSim: RDWaveform
Base.:(==)(wf::RDWaveform, nwf::RDWaveform) = 
    reduce((x,field) -> x && isequal(getfield(wf, field), getfield(nwf, field)), fieldnames(RDWaveform), init = typeof(wf) == typeof(nwf))


@testset "LegendGeSim: Simulate capacitances of test detectors" begin

    testdata_path = joinpath(legend_test_data_path(), "data", "legend", "metadata", "hardware", "detectors", "germanium", "diodes")
    test_dict = Dict{String, typeof(1.0u"pF")}(
        "B99000A.json" => -7.13u"pF",
        "V99000A.json" => -3.55u"pF"
    )

    for (filename, capacitance) in test_dict
        detector_metadata_filename = joinpath(testdata_path, filename)
        sim_settings_ssd_filename = joinpath(dirname(dirname(pathof(LegendGeSim))), "examples/configs/detector_study_ssd.json")
        
        @testset "Field Simulation of $filename (SSD)" begin
            sim = LegendGeSim.simulate_fields(detector_metadata_filename, sim_settings_ssd_filename)
            C = LegendGeSim.capacitance_matrix(sim)

            @testset "Capacitances" begin   
                @test isapprox(C[1,2], capacitance, atol = 0.2u"pF")
            end
        end
    end 
end


@testset "LegendGeSim: sim -> raw simulation chain" begin

    ldsim_path = joinpath(legend_test_data_path(), "data", "ldsim")
    det_metadata = joinpath(ldsim_path, "invcoax-metadata.json")
    path_to_pet_file = joinpath(ldsim_path, "single-invcoax-th228-geant4.csv")

    environment_settings = Dict(
        "crystal_temperature_in_K" => 77,
        "medium" => "vacuum" 
    )

    simulation_settings = Dict(
        "method" => "ssd",
        "cached_name" => "test",
    )
    

    pss_table, pss_truth = LegendGeSim.simulate_pulses(det_metadata, path_to_pet_file, environment_settings, simulation_settings) #; n_waveforms=10

    @testset "pet -> pss" begin 
        @test pss_table isa LegendGeSim.Table
        @test pss_truth isa LegendGeSim.Table
    end

    setup_settings = Dict(
        "preamp" => Dict(
            "type" => "generic",
            "t_decay" => 50,
            "t_rise" => 15,
            "max_e" => 10000,
            "offset" => 2000
        ),
        "fadc" => Dict(
            "type" => "generic",
            "sampling_interval" => 16
        ),
        "trigger" => Dict(
            "type" => "trapezoidal",        
            "window_lengths" => [250,250,250],
            "threshold" => 9
        ),
        "daq" => Dict(
            "type" => "generic",
            "nsamples" => 3750,
            "baseline_length" => 1875
        )
    )

    pss_name = "cache/" * LegendGeSim.filename(path_to_pet_file) * "_pss.hdf5"
    h5open(pss_name, "w") do f
        LegendHDF5IO.writedata(f, "pss/pss", pss_table)
        LegendHDF5IO.writedata(f, "pss/truth", pss_truth)
    end

    @testset "pss I/O using LegendHDF5IO" begin
        @test isfile(pss_name)
        h = LHDataStore(pss_name)
        @test haskey(h, "pss")
        close(h)
    end

    raw_table = LegendGeSim.pss_to_raw(pss_name, setup_settings)

    @testset "pss -> raw" begin
        @test raw_table isa LegendGeSim.Table
        @test length(raw_table) == 240
    end

    # Skip writing the pss to disk and check that the result is still the same
    all_settings = Dict(
        "environment" => environment_settings,
        "simulation" => simulation_settings,
        "setup" => setup_settings
    )

    new_raw_table = LegendGeSim.simulate_raw(det_metadata, path_to_pet_file, all_settings; n_waveforms=5)

    @testset "pet -> raw" begin
        @test new_raw_table[1:5] == raw_table[1:5]
    end

    # Delete the cached files
    rm("cache", force=true, recursive=true)
end
