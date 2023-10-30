# This file is a part of LegendGeSim.jl, licensed under the MIT License (MIT).

to_SSD_units(::Type{T}, x, unit) where {T} = T(SolidStateDetectors.to_internal_units(x*unit)) 

# if just want geometry
function LEGEND_SolidStateDetector(metapath::AbstractString)
    LEGEND_SolidStateDetector(propdict(metapath))
end

function LEGEND_SolidStateDetector(meta::PropDict)
    SolidStateDetector{Float32}(LegendData, meta)
end

function LEGEND_SolidStateDetector(::Type{T}, meta::PropDict, env::Environment, impurity_path::AbstractString = "") where {T}
    # ------------------------------------------------------------------------------
    # temporary DL quickfix through Environment - super ugly, needs to be changed 
    
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

    # add field geometry.dl_thickness_in_mm to meta dict - LegendDataManagement will look for that to set DL 
    meta.geometry[:dl_thickness_in_mm] = dl_thickness_in_mm

    # ------------------------------------------------------------------------------
    # geometry

    detector = SolidStateDetector{T}(LegendData, meta)

    # ------------------------------------------------------------------------------
    # operating voltage

    # sometimes recV is null -> in this case, complain and ask to provide (default env.opV is 0 if none given)
    operation_voltage = 
    if env.operating_voltage == 0
        @warn "You did not provide operating voltage in environment settings -> taking recommended voltage from metadata"
        recV = meta.characterization.l200_site.recommended_voltage_in_V
        if isnothing(recV)
            error("Metadata does not provide recommended voltage (null). Please provide operating voltage in settings")
        end
        T(recV)
    else
        T(env.operating_voltage)
    end

    @info "Simulating at $(operation_voltage)V"

    detector = SolidStateDetector(detector, contact_id = 2, contact_potential = operation_voltage)

    # ------------------------------------------------------------------------------
    # temperature

    temperature = T(env.crystal_temperature)
    # "old" semiconductor 
    sc = detector.semiconductor
    # new semiconductor with our temperature
    semiconductor = SolidStateDetectors.Semiconductor(temperature, sc.material, sc.impurity_density_model, sc.charge_drift_model, sc.geometry)
    # new detector with new semiconductor
    detector = SolidStateDetector{T}(detector.name, semiconductor, detector.contacts, detector.passives, detector.virtual_drift_volumes)

    # ------------------------------------------------------------------------------
    # impurity profile

    # if no impurity path given, simulate with constant impurity density
    imp_density_model = 
        if impurity_path == ""
            # Impurity Model: Information are stored in `meta.production.impcc`
            # For now: Constant impurity density: 
            #   n-type: positive impurity density
            #   p-type: negative impurity density
            # Assume p-type
            constant_impurity_density = ustrip(uconvert(u"m^-3", T(-9*1e9) * u"cm^-3"))
            SolidStateDetectors.CylindricalImpurityDensity{T}(
                (0, 0, constant_impurity_density), # offsets
                (0, 0, 0)                          # linear slopes
            )
        else
            # currently using David Radford's empirical function
            # quadratic implemented as well but 1) David's is better, and 2) haven't implemented keywords to let user choose
            impurity_density_model(meta, impurity_path, SSDSimulator())
        end

    detector = SolidStateDetector(detector, imp_density_model)
end