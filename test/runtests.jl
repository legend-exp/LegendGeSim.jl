# This file is a part of LegendGeSim.jl, licensed under the MIT License (MIT).

using Test

using LegendGeSim
using LegendHDF5IO, HDF5
using LegendTestData
using Unitful


@testset "LegendGeSim: Simulate capacitances of test detectors" begin
    # note: this part also tests dict settings input as opposed to json file
    # as well as impurity input
    germanium_testdata_path = joinpath(legend_test_data_path(), "data", "legend", "metadata", "hardware", "detectors", "germanium")
    test_dict = Dict{String, typeof(1.0u"pF")}(
        "B99000A.json" => -5.20u"pF",
        "V99000A.json" => -3.55u"pF"
    )

    environment_settings = Dict(
        "crystal_temperature_in_K" => 77,
        "medium" => "vacuum" 
    )

    simulation_settings = Dict(
        "crystal_metadata_path" => joinpath(germanium_testdata_path, "crystals")
    )


    for (filename, capacitance) in test_dict
        detector_metadata_filename = joinpath(germanium_testdata_path, "diodes", filename)

        for method in ["ssd", "siggen"]
            simulation_settings["method"] = method

            @testset "Field Simulation of $filename ($method)" begin
                sim = LegendGeSim.simulate_fields(detector_metadata_filename, environment_settings, simulation_settings)
                C = LegendGeSim.capacitance_matrix(sim)

                @testset "Capacitances" begin   
                    @test isapprox(C[1,2], capacitance, atol = 0.4u"pF")
                end
            end
        end
    end 
end


@testset "LegendGeSim: sim -> raw simulation chain" begin
    # note: this part also tests json settings input as opposed to dict
    # as well as simulation with no impurity input i.e. dummy constant impurity
    ldsim_path = joinpath(legend_test_data_path(), "data", "ldsim")
    det_metadata = joinpath(ldsim_path, "invcoax-metadata.json")
    path_to_pet_file = joinpath(ldsim_path, "single-invcoax-th228-geant4.csv")
    sim_settings_ssd_filename = joinpath(dirname(dirname(pathof(LegendGeSim))), "examples/configs/detector_study_ssd.json")

    pss_table, pss_truth = LegendGeSim.simulate_pulses(det_metadata, path_to_pet_file, sim_settings_ssd_filename) #; n_waveforms=10

    @testset "pet -> pss" begin 
        @test pss_table isa LegendGeSim.Table
        @test pss_truth isa LegendGeSim.Table
    end

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
    
    # ToDo: check eV vs ADC parameters: 1) either should work, 2) error when both provided
    # ToDo: if noise model is "none" but preamp sigma is given, will simulate noise anyway! need to fix
    setup_settings = Dict(
        "preamp" => Dict(
            "type" => "generic",
            "t_decay_in_us" => 50,
            "t_rise_in_ns" => 100,
            "gain_ADC_eV" => 0.138,
            "offset_in_ADC" => 12000
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

    raw_table = LegendGeSim.pss_to_raw(pss_name, setup_settings)

    @testset "pss -> raw (without noise)" begin
        @test raw_table isa LegendGeSim.Table
        @test length(raw_table) == 240
        @test sum(!iszero, raw_table.energy) == 238 # two do not trigger
    end
    
    # note: this is a dummy test of electronics noise
    # since no noise model was provided for pet -> pss, there is no fano noise at that stage
    # so this is not truly the NoiseFromSim simulation 
    # Note: there is no way to check on user level whether the previous pre-simulated pet->pss had fano noise -> baseline std?
    # ToDo: noise from data [still very WIP so too early to test]
    setup_settings_noise = deepcopy(setup_settings)
    setup_settings_noise["preamp"]["noise_sigma_in_keV"] = 2
    noise_settings = Dict( "type" => "sim" )
    
    raw_table_noise = LegendGeSim.pss_to_raw(pss_name, setup_settings_noise, noise_settings)
    
    @testset "pss -> raw (with noise)" begin
        @test raw_table isa LegendGeSim.Table
        @test length(raw_table) == 240
        @test sum(!iszero, raw_table.energy) == 238 # two do not trigger
    end


    # Skip writing the pss to disk and check that the result is still the same
    all_settings = Dict(
        "environment" => Dict(
            "crystal_temperature_in_K" => 77,
            "medium" => "vacuum" 
        ),
        "simulation" =>  Dict(
            "method" => "SSD"
            # note: there shouldn't be crystal metadata here, as the json used in pet->pss doesn't have it 
            # note: to not have bugs because two different sources of settings, maybe avoid this
        ),
        "setup" => setup_settings
    )

    new_raw_table = LegendGeSim.simulate_raw(det_metadata, path_to_pet_file, all_settings; n_waveforms=5)

    @testset "pet -> raw" begin
        @test new_raw_table[1:5] == raw_table[1:5]
    end

    # Delete the cached files
    rm("cache", force=true, recursive=true)
end
