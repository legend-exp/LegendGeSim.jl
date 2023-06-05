"""
    siggen_config(meta, env, siggen_sim, config_name)

PropDict, Environment, SiggenSimulator, String -> String, String    

Construct Siggen/Fieldgen config based on
    geometry as given in LEGEND metadata <meta>,
    environment settings given in <env>,
    and fieldgen settings contained in <siggen_sim>.

The resulting siggen config will be cached based on given <config_name>
The function returns the path to the cached siggen config file, and
    the relative path to the fieldgen WP file (to check if already exists later)
"""
function siggen_config(meta::PropDict, env::Environment, simulator::SiggenSimulator; overwrite::Bool=false)
    # ----------------------------------------------------------------------------
    # set up & check if already present in cache
    # ----------------------------------------------------------------------------

    # if no cached name was given by user, make a temporaty name - needed for fieldgen to read from file
    # also if no cached name was given means definitely doing from scratch
    if simulator.cached_name == ""
        simulator.cached_name = "tmp"
        overwrite = true
    end

    # cached_name = simulator.cached_name == "" ? "tmp" : simulator.cached_name    
    # append detector name
    cached_name = meta.name * "_" * simulator.cached_name  
    
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

    # why is this needed? 
    for name in ["field_name", "wp_name"]
        path = joinpath("cache", fieldgen_names[name])
        if !isfile(path)
            touch(path)
        end
    end

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


"""
    meta2siggen(meta, env)

PropDict, Environment -> Dict    

Construct siggen type config in a form of Dict based on
    geometry as given in LEGEND metadata <meta>
    and environment settings given in <env>.
"""
function meta2siggen(meta::PropDict, env::Environment)
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
        # depth of full-charge-collection boundary for Li contact -> was removed from geometry because is not 
        # use manufacturer DL? in the future when our FCCD is in metadata, use that?
        # "Li_thickness" => 0,
        "Li_thickness" => meta.characterization.manufacturer.dl_thickness_in_mm,

        # configuration for mjd_fieldgen (calculates electric fields & weighing potentials)
        # detector bias for fieldgen, in Volts
        # ToDo: sometime there is null in metadata -> throw error and ask for opV input
        "xtal_HV" => env.operating_voltage > 0 ? env.operating_voltage : meta.characterization.l200_site.recommended_voltage_in_V,

        # configuration for signal calculation 
        # crystal temperature in Kelvin
        "xtal_temp" => env.crystal_temperature
    )

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