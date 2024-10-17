module LegendGeSimExt
#using MJDSigGen: SiggenSimulator
#@info "Error here"
using LegendGeSim
#@info "Error"
using PropDicts
#@info "Error_prop"
using LegendGeSim
#@info "Error_Lgsim"
using LegendGeSim: Environment
#@info "le_env"
#using MJDSigGen: SiggenSimulator
using MJDSigGen
using MJDSigGen: SigGenSetup
#@info "Error_mjd"
#using LegendGeSim: PSSimulator
using LegendGeSim: SiggenSimulator
using ArgCheck
using ArraysOfArrays
using CurveFit
using DelimitedFiles
using Distributions
using DSP 
using ElasticArrays
using HDF5
using IntervalSets
using JSON
using LegendDataManagement
using LegendDataTypes
using LegendHDF5IO
using LegendTextIO
using LinearAlgebra
using Parameters
using Polynomials
using PropDicts
using RadiationDetectorDSP 
using RadiationDetectorSignals
using RadiationSpectra
using Random
using Random123
using RecipesBase
using Requires
using SignalAnalysis
using SolidStateDetectors
using StaticArrays
using Statistics
using StatsBase
using StructArrays
using Tables
using TypedTables
using Unitful

import LsqFit

using SolidStateDetectors: ConstructiveSolidGeometry as CSG

#abstract type PSSimulator end


#include("legend_detector_to_siggen.jl")
#@info "Here_then"
#include("mjdsiggen_utils.jl")
#@info "Error_after_files"



####################################
### FieldGen
####################################
@info "27"
"""
    SiggeSimulator(sim_conf::PropDict)

    LegendGeSimConfig -> SiggenSimulator

Construct SiggeSimulator instance based on simulation
    configuration given in <sim_conf>.
"""

function SiggenSimulator(simulation_settings::PropDict)
    # set defaults
    inputs = Dict{String, Any}()
    for k in (:fieldgen_config, :drift_vel, :impurity_profile, :crystal_metadata_path)
        inputs[string(k)] = haskey(simulation_settings, k) ? simulation_settings[k] : ""
    end    

    # check if provided paths exist
    for (param, path) in inputs
        if (path != "") && !(isfile(path) || isdir(path))
        throw(ArgumentError("The input for $(param) that you provided does not exist: $(path)"))
        end
    end

    inputs["offset_in_mm"] = haskey(simulation_settings, :offset_in_mm) ? simulation_settings.offset_in_mm : -1

    # throw error if both files are provided
    if inputs["crystal_metadata_path"] != "" && inputs["impurity_profile"] != ""
        throw(ArgumentError("You provided both an .spe/.dat file and path to crystal metadata! Provide one or the other as impurity input for siggen."))
    # warn if none are provided - will use constant impurity density
    elseif inputs["crystal_metadata_path"] == "" && inputs["impurity_profile"] == ""
        @warn "No impurity inputs given. Simulation with dummy constant impurity density."
    # if we get here, one is provided
    else
        impinput = inputs["impurity_profile"] * inputs["crystal_metadata_path"]
        @info "Impurity profile information based on $(impinput)"
    end

    # check that offset is provided if .spe/.dat file is
    if inputs["impurity_profile"] != "" && inputs["offset_in_mm"] == -1
        throw(ArgumentError("Provide offset_in_mm of this detector corresponding to impurity file $(inputs["impurity_profile"])!"))
    end


    SiggenSimulator( 
        inputs["fieldgen_config"],
        inputs["drift_vel"],
        inputs["impurity_profile"],
        inputs["offset_in_mm"],
        inputs["crystal_metadata_path"],
        simulation_settings.cached_name
    )
end

"""
    simulate_detector(det_meta, env, config_name, ps_simulator)
    
PropDict, Environment, AbstractString, Siggen -> SigGenSetup       

Construct a fieldgen/siggen configuration file
    according to geometry as given in LEGEND metadata <det_meta>
    and environment settings given in <env>.
Look up fieldgen generated electric field and weighting potential files
    based on given <config_name>, and if not found, call fieldgen.
"""

function LegendGeSim.simulate_detector(det_meta::PropDict, env::Environment, simulator::SiggenSimulator;
    overwrite::Bool = false)

# if no cached name given, force simulation from scratch (or else might read past "tmp" file)
if simulator.cached_name == ""
    overwrite = true
end

# returns the name of the resulting siggen config file
# and the name of (already or to be) generated weighting potential file
# ToDo: don't create if already present already at this stage -> check in siggen_config() if name exists and just return name
siggen_config_name, fieldgen_wp_name = siggen_config(det_meta, env, simulator; overwrite)
fieldgen_wp_name = joinpath("cache", fieldgen_wp_name)

# if the WP file with such a name does not exist yet or want to overwrite it...
if !isfile(fieldgen_wp_name) || overwrite
    imp_filename, offset =
        if simulator.crystal_metadata_path != ""
            # create impurity input on the fly
            LegendGeSim.impurity_density_model(det_meta, simulator.crystal_metadata_path, simulator)
        else
            # user provided .dat file and corresponding offset
            # (or if nothing is provided this will be "" and -1)
            simulator.impurity_profile, simulator.offset_in_mm
        end
    #...call fieldgen -> will write the wp file
    @info "_|~|_|~|_|~|_ Fieldgen simulation"
    fieldgen(siggen_config_name; impurity_profile=imp_filename, offset_mm=offset)
    @info "_|~|_|~|_|~|_ Fieldgen simulation complete"
else
    #...do nothing, siggen will later read the files based on the generated conifg
    @info "Cached simulation found. Reading cached fields from $fieldgen_wp_name"
end
@info  "73"
# a SigGenSetup object
SigGenSetup(siggen_config_name)
end

## ---------- fieldgen

"""
Convert fit parameters from crystal metadata units (1e9 e/cm^3 VS mm) into Siggen units (1e10 e/cm^3 VS mm)
and create a .dat file (unformatted stream of Float32) to be later used when calling fieldgen().
"""
function LegendGeSim.impurity_density_model(meta::PropDict, crystal_metadata::PropDict, fit_param::Vector{Float64}, ::SiggenSimulator)   

    a,b,c,tau = fit_param 

    T = Float32
    dist = T.(crystal_metadata.impurity_measurements.distance_from_seed_end_mm)
    # sometimes might not have measurements at seed end? but last measurement can be taken as crystal length
    L = dist[end] #mm

    # points in crystal axis from 0 to L with step of 1 mm
    # always put 200 points to make all .dat files same length
    imp = zeros(Float32, 200) # impurities in 10^10 cm^-3

    # length is L+1 including both ends
    Npoints = Int(ceil(L))+1 # step of 1 mm
    zpoints = LinRange(0, L, Npoints) 
    for i in 1:Npoints
        # divide by 10 to convert from 1e9 e/cm^3 to 1e10 e/cm^3
        # invert the order of the values because with siggen z=0 is tail (dirtiest part)
        # i.e. z = L in our system (Mirion measurements)
        imp[Npoints-i+1] = -(a + b * zpoints[i] + c * exp((zpoints[i] - L)/tau)) / 10.
    end

    crystal_name = first(meta.name, length(meta.name)-1)
    imp_filename = joinpath("cache", crystal_name*".dat")
	open(imp_filename, "w") do io
		write(io, imp)
	end    
    @info "Impurity profile for siggen saved in $(imp_filename)"

    # return name of impurity profile file and offset
    # offset from seed end
    det_z0 = crystal_metadata.slices[Symbol(meta.production.slice)].detector_offset_in_mm    

    # siggen offset is from tail i.e. at L offset = 0
    imp_filename, L - det_z0
end

@info "123"
"""
Simulation method: siggen
"""



function LegendGeSim.simulate_waveforms(stp_events::Table, detector::SigGenSetup, ::SiggenSimulator)
    T = Float32 # This should be somehow defined and be passed properly
    println("127")
    @info("~.~.~ Siggen")
    @info "129"

    nevents = length(stp_events)
    wf_array = Array{RDWaveform}(undef, nevents)
    for i in 1:nevents
        # println("$i/$nevents")
        if(i % 100  == 0) println("$i/$nevents") end
        signal::Vector{Float32} = simulate_signal(stp_events[i].pos, stp_events[i].edep, detector)
        time = range(T(0)u"ns", step = detector.step_time_out*u"ns", length = length(signal))

        # push!(waveforms, signal)
        wf_array[i] = RDWaveform(time, signal)
    end

    ## Oliver said this should be the more efficient way of doing it
    ## as opposed to an array of separate RDWaveform instances
    # Δt = setup.step_time_out*u"ns"
    # t_end = length(waveforms[1]) * Δt
    # times = fill(0u"ns" : Δt : t_end, length(waveforms))
    # ArrayOfRDWaveforms((times, VectorOfSimilarVectors(waveforms)))

    waveforms = ArrayOfRDWaveforms(wf_array)
    # why am I doing this?
    waveforms = ArrayOfRDWaveforms((waveforms.time, VectorOfSimilarVectors(waveforms.signal)))    

    pss_table = Table(
        channel = [1 for idx in 1:length(waveforms)], # lists of ADCs that triggered, 1 for HADES all the time
        ievt = stp_events.evtno,
        waveform = waveforms
    )

    # using SSD we get mc truth from SSD output
    # with siggen, let's just return the input
    # maybe mc truth will not be propagated to (sim)raw anyway
    pss_truth = Table(
        detno = stp_events.detno,
        edep = stp_events.edep,
        ievt = stp_events.evtno,
        pos = stp_events.pos,
        thit = stp_events.thit        
    )

    pss_table, pss_truth
end


#@info "252"

"""
    simulate_wf(pos, edep, siggen_setup)

AbstractVector, AbstractVector, SigGenSetup -> Vector{Float32}

Simulate a signal from energy depositions <edep> with corresponding positions <pos>
    based on a given <siggen_setup>    
"""

function simulate_signal(pos::AbstractVector, edep::AbstractVector, siggen_setup::SigGenSetup)
    # this is not so nice, think about it
    signal = zeros(Float32, siggen_setup.ntsteps_out)

    # round to siggen crystal grid precision
    # a = string(siggen_setup.xtal_grid) # crystal grid precision e.g. 0.1 (mm)
    # ndig = length(a[findfirst('.', a)+1:end]) # number of digits after .

    for i in 1:length(pos)
        # pos_rounded = round.(ustrip.(pos[i]), digits=ndig) # strip of units and round to siggen precision
        # in SSD the output is always in eV -> convert to eV
        signal = signal .+ MJDSigGen.get_signal!(siggen_setup, Tuple(ustrip.(pos[i]))) * ustrip(uconvert(u"eV", edep[i]))
        # signal = signal .+ MJDSigGen.get_signal!(siggen_setup, Tuple(pos_rounded)) * ustrip(uconvert(u"eV", edep[i]))

        # somehow this doesn't work
        # MJDSigGen.get_signal!(signal, siggen_setup, Tuple(ustrip.(pos[i]))) * ustrip(uconvert(u"eV", edep[i]))
    end

    signal
end

#Calculating capacitance matrix using Siggen 
function LegendGeSim.capacitance_matrix(sim::SigGenSetup) 
    c = sim.capacitance * u"pF"
    [c -c; -c c]
end
#export capacitance_matrix
#@info "288"
#legend_detector_to_siggen
function siggen_config(meta::PropDict, env::Environment, simulator::SiggenSimulator; overwrite::Bool=false)
    # ----------------------------------------------------------------------------
    # set up & check if already present in cache
    # ----------------------------------------------------------------------------

    # if no cached name was given by user, make a temporaty name - needed for fieldgen to read from file
    # also if no cached name was given means definitely doing from scratch
    cached_name = simulator.cached_name
    if simulator.cached_name == ""
        cached_name = "tmp"
        overwrite = true
    end

    # append detector name
    cached_name = meta.name * "_" * cached_name  
    
    # filenames for fieldgen input/output 
    # ! no need to add cache/ folder because siggen does it by itself...
    fieldgen_names = Dict(
        "field_name" => cached_name * "_fieldgen_EV.dat",
        "wp_name" => cached_name * "_fieldgen_WP.dat"
    )

    # siggen config filename
    sig_conf_name = joinpath("cache", cached_name * "_siggen_config.txt")

    # if field exists, no need to construct config, just read cached simulation -> return
    # if config exists, no need to construct again -> return -> maybe just the cached field? not sure
    # of course unless asked to overwrite, then construct new config
    if ( isfile(joinpath("cache", fieldgen_names["wp_name"])) || isfile(sig_conf_name) ) && !overwrite
        return sig_conf_name, fieldgen_names["wp_name"]
    end

    @info "Constructing Fieldgen/Siggen configuration file"

    if !ispath("cache") mkpath("cache") end


    # ----------------------------------------------------------------------------
    # construct config
    # ----------------------------------------------------------------------------

    println("...detector geometry")
    
    # construct siggen geometry based on LEGEND metadata JSON + environment settings
    siggen_geom = meta2siggen(meta, env)

    # collect geometry lines for siggen config file 
    detector_lines = Vector{String}()
    push!(detector_lines, "# detector related parameters")
    for (param, value) in siggen_geom
        push!(detector_lines, "$param   $value")
    end

    println("...fieldgen configuration")

    # ----------------------------------------------------------------------------
    # drift velocity, field, wp
    # ----------------------------------------------------------------------------

    # drift velocity input (path to file given by user)
    drift_file_path = simulator.drift_vel
    # use default if not given by user
    if simulator.drift_vel == ""
        @warn "No drift velocity input given; using default"
        drift_file_path = joinpath(dirname(@__DIR__), "examples", "configs", "drift_vel_tcorr.tab")      
    end

    # copy to cache
    cache_drift_file = cached_name * "_" * basename(drift_file_path)
    cp(drift_file_path, joinpath("cache", cache_drift_file), force=true) # later think what to do here

    # add to field fieldgen_names
    fieldgen_names["drift_name"] = cache_drift_file

    ## commenting this out: produced a bug
    ## (empty file created, next time running empty file read and no field at point error)
    # for name in ["field_name", "wp_name"]
    #     path = joinpath("cache", fieldgen_names[name])
    #     if !isfile(path)
    #         touch(path)
    #     end
    # end

    # add corresponding lines to existing detector geometry lines
    push!(detector_lines, "")
    push!(detector_lines, "# fieldgen filenames")
    for (param, value) in fieldgen_names
        push!(detector_lines, "$param   $value")
    end

    # ----------------------------------------------------------------------------
    # extra settings
    # ----------------------------------------------------------------------------

    # Note: this is similar but not identical to the part above because of the cache folder difference
    # -> loop? copy-paste is not nice...

    # read lines from user input for fieldgen
    config_file_path = simulator.fieldgen_config
    if simulator.fieldgen_config == ""
        @warn "No extra fieldgen settings given; using default"
        config_file_path = joinpath(dirname(@__DIR__), "examples", "configs", "fieldgen_settings.txt")      
    end
    # copy to cache -> maybe not needed, are appended in main settings?
    cache_config_file = joinpath("cache", cached_name * "_" * basename(config_file_path))
    cp(config_file_path, cache_config_file, force=true) # later think what to do here

    fieldgen_lines = readlines(open(cache_config_file))
    fieldgen_lines = replace.(fieldgen_lines, "\t" => "   ")

    # unite constructed detector geometry and fieldgen settings
    total_lines = vcat(detector_lines, [""], fieldgen_lines)

    # ----------------------------------------------------------------------------
    # final siggen config file
    # ----------------------------------------------------------------------------

    # write a siggen/fieldgen config file 

    writedlm(sig_conf_name, total_lines)
    println("...fieldgen/siggen config written to $sig_conf_name")

    # return resulting fieldgen/siggen config name 
    # (currently kinda redundant but that's because we have no idea what's gonna happen later)
    sig_conf_name, fieldgen_names["wp_name"]
end
#@info "417"

"""
    meta2siggen(meta, env)

PropDict, Environment -> Dict    

Construct siggen type config in a form of Dict based on
    geometry as given in LEGEND metadata <meta>
    and environment settings given in <env>.
"""
function meta2siggen(meta::PropDict, env::Environment)
    # !!!! this is copy paste from legend_detector_to_SSD !!
    # ToDo: rewrite to avoid copy-paste 

    # -----------------------------------------------------------
    # DL thickness (quickfix)
        
    dl_thickness_in_mm =
        if env.dl == "vendor"
            dl_vendor = meta.characterization.manufacturer.dl_thickness_in_mm
            if dl_vendor isa PropDicts.MissingProperty
                throw(ArgumentError("No dead layer measurement provided by vendor! Please provide value or skip (default 0)."))
            end
            dl_vendor
        else
            env.dl 
        end

    @info "DL = $(dl_thickness_in_mm)mm"

    # -----------------------------------------------------------
    # operating voltage

    # sometimes recV is null -> in this case, complain and ask to provide (default env.opV is 0 if none given)
    operation_voltage = 
    if env.operating_voltage == 0
        @warn "You did not provide operating voltage in environment settings -> taking recommended voltage from metadata"
        recV = meta.characterization.l200_site.recommended_voltage_in_V
        if isnothing(recV)
            error("Metadata does not provide recommended voltage (null). Please provide operating voltage in settings")
        end
        recV
    else
        env.operating_voltage
    end

    @info "Simulating at $(operation_voltage)V"    

    # -----------------------------------------------------------
    # parameters common for all geometries

    siggen_geometry = Dict(
        # -- crystal
        "xtal_length" => meta.geometry.height_in_mm,
        "xtal_radius" => meta.geometry.radius_in_mm,
        # -- p+ contact
        "pc_length" => meta.geometry.pp_contact.depth_in_mm,
        "pc_radius" => meta.geometry.pp_contact.radius_in_mm,
        # -- tapering

        # comment from David Radford:

        # The bottom taper is always at 45 degrees, but the top outer taper is not.
        # It's width/angle is defined using either the outer_taper_width parameter or the taper_angle parameter.
        # Likewise, the inner taper width/angle is defined using either the inner_taper_width parameter or the taper_angle parameter.
        # The inner_taper_length can be anything from zero to the hole_length

        # size of 45-degree taper at bottom of ORTEC-type crystal (for r=z)
        # -> for now, ignore bottom taper if not 45
        # -> works for 1) the ones with actual 45 taper, 2) the ones bulletized (saved as 45 deg taper with bullet radius)
        # few that have proper non-45 taper - not possible. ToDo: put some warning
        "bottom_taper_length" => meta.geometry.taper.bottom.angle_in_deg == 45 ? meta.geometry.taper.bottom.height_in_mm : 0,
        # current metadata format: always angle
        "bottom_taper_width" => meta.geometry.taper.bottom.angle_in_deg == 45 ? meta.geometry.taper.bottom.height_in_mm * tan(deg2rad(meta.geometry.taper.bottom.angle_in_deg)) : 0,
        # z-length of outside taper for inverted-coax style -> Mariia: you mean top taper here?
        "outer_taper_length" => meta.geometry.taper.top.height_in_mm,
        "outer_taper_width" => meta.geometry.taper.top.height_in_mm * tan(deg2rad(meta.geometry.taper.top.angle_in_deg)),
        # taper angle in degrees, for inner or outer taper -> how? why are they the same? how about borehole? not using this parameter...
        # "taper_angle" => meta.geometry.taper.top.outer.angle_in_deg,
        # "taper_angle" => meta.geometry.taper.top.inner.angle_in_deg,
        # with this setting, it somehow takes taper angle: 45.000038 anyway -> force to zero?
        # "taper_angle" => 0,
        # depth of full-charge-collection boundary for Li contact
        "Li_thickness" => dl_thickness_in_mm,

        # configuration for mjd_fieldgen (calculates electric fields & weighing potentials)
        # detector bias for fieldgen, in Volts
        "xtal_HV" => operation_voltage,

        # configuration for signal calculation 
        # crystal temperature in Kelvin
        "xtal_temp" => env.crystal_temperature
    )

    # -----------------------------------------------------------
    # non common parameters

    # -- groove (for non PPCs)
    if haskey(meta.geometry, :groove)
        siggen_geometry["wrap_around_radius"] = meta.geometry.groove.radius_in_mm.outer
        siggen_geometry["ditch_depth"] = meta.geometry.groove.depth_in_mm
        siggen_geometry["ditch_thickness"] = meta.geometry.groove.radius_in_mm.outer - meta.geometry.groove.radius_in_mm.inner
    end

    # -- borehole (for ICPC and Coax)
    if haskey(meta.geometry, :borehole)
        # length of hole, for inverted-coax style
        siggen_geometry["hole_length"] = meta.geometry.borehole.depth_in_mm
        # radius of hole, for inverted-coax style
        siggen_geometry["hole_radius"] = meta.geometry.borehole.radius_in_mm
        # z-length of inside (hole) taper for inverted-coax style
        siggen_geometry["inner_taper_length"] = meta.geometry.taper.borehole.height_in_mm
        siggen_geometry["inner_taper_width"] = meta.geometry.taper.borehole.height_in_mm * tan(deg2rad(meta.geometry.taper.borehole.angle_in_deg))
    end

    siggen_geometry
end

#mjdsiggen_utils.jl
#@info "537"

function SolidStateDetectors.ElectricPotential(sim::SigGenSetup)
    E_pot, W_pot, E_abs, E_r, E_z = LegendGeSim.MJDSigGen.read_fields(sim);
    T = eltype(E_pot)
    r_axis = (0:(sim.rlen - 1)) * sim.xtal_grid
    z_axis = (0:(sim.zlen - 1)) * sim.xtal_grid
    ax1 = SolidStateDetectors.DiscreteAxis{T, :r0, :fixed, ClosedInterval{T}}(r_axis[1]..r_axis[end], r_axis)
    ax2 = SolidStateDetectors.DiscreteAxis{T, :reflecting, :reflecting, ClosedInterval{T}}(zero(T)..zero(T), zeros(T, 1))
    ax3 = SolidStateDetectors.DiscreteAxis{T, :reflecting, :reflecting, ClosedInterval{T}}(z_axis[1]..z_axis[end], z_axis)
    axs = (ax1, ax2, ax3)
    grid = Grid{T, 3, SolidStateDetectors.Cylindrical, typeof(axs)}( axs )
    data = Array{T, 3}(undef, length(ax1), length(ax2), length(ax3));
    data[:] = E_pot'[:];
    ElectricPotential(data, grid)
end

function SolidStateDetectors.WeightingPotential(sim::SigGenSetup)
    E_pot, W_pot, E_abs, E_r, E_z = LegendGeSim.MJDSigGen.read_fields(sim);
    T = eltype(W_pot)
    r_axis = (0:(sim.rlen - 1)) * sim.xtal_grid
    z_axis = (0:(sim.zlen - 1)) * sim.xtal_grid
    ax1 = SolidStateDetectors.DiscreteAxis{T, :r0, :fixed, ClosedInterval{T}}(r_axis[1]..r_axis[end], r_axis)
    ax2 = SolidStateDetectors.DiscreteAxis{T, :reflecting, :reflecting, ClosedInterval{T}}(zero(T)..zero(T), zeros(T, 1))
    ax3 = SolidStateDetectors.DiscreteAxis{T, :reflecting, :reflecting, ClosedInterval{T}}(z_axis[1]..z_axis[end], z_axis)
    axs = (ax1, ax2, ax3)
    grid = Grid{T, 3, SolidStateDetectors.Cylindrical, typeof(axs)}( axs )
    data = Array{T, 3}(undef, length(ax1), length(ax2), length(ax3));
    data[:] = W_pot'[:];
    WeightingPotential(data, grid)
end

end #module LegendGeSimExt
