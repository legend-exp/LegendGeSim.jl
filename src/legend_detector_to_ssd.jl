# This file is a part of LegendGeSim.jl, licensed under the MIT License (MIT).

to_SSD_units(::Type{T}, x, unit) where {T} = T(SolidStateDetectors.to_internal_units(x*unit)) 

# if just want geometry
function LEGEND_SolidStateDetector(metapath::AbstractString)
    LEGEND_SolidStateDetector(Float32, propdict(metapath), Environment("",0,0))
end

function LEGEND_SolidStateDetector(::Type{T}, meta::PropDict, env::Environment, impurity_path::AbstractString = "") where {T}
        # Not all possible configurations are yet implemented!
        # https://github.com/legend-exp/legend-metadata/blob/main/hardware/detectors/detector-metadata_1.pdf
        # https://github.com/legend-exp/legend-metadata/blob/main/hardware/detectors/detector-metadata_2.pdf
        # https://github.com/legend-exp/legend-metadata/blob/main/hardware/detectors/detector-metadata_3.pdf
        # https://github.com/legend-exp/legend-metadata/blob/main/hardware/detectors/detector-metadata_4.pdf
        # https://github.com/legend-exp/legend-metadata/blob/main/hardware/detectors/detector-metadata_5.pdf
        # https://github.com/legend-exp/legend-metadata/blob/main/hardware/detectors/detector-metadata_6.pdf
        # https://github.com/legend-exp/legend-metadata/blob/main/hardware/detectors/detector-metadata_7.pdf

        gap = to_SSD_units(T, 1, u"mm")

        # dl_thickness_in_mm = :dl_thickness_in_mm in keys(meta.geometry) ? meta.geometry.dl_thickness_in_mm : 0
        # ToDo: read manufacturer DL value. DL was removed from geometry
        dl_thickness_in_mm = 0
        li_thickness =  to_SSD_units(T, dl_thickness_in_mm, u"mm")

        crystal_radius = to_SSD_units(T, meta.geometry.radius_in_mm, u"mm")
        crystal_height = to_SSD_units(T, meta.geometry.height_in_mm, u"mm")

        # main crystal
        semiconductor_geometry = CSG.Cone{T}(CSG.ClosedPrimitive; 
            r = crystal_radius, 
            hZ = crystal_height / 2, 
            origin = CartesianPoint{T}(0, 0, crystal_height / 2)
        )
            
        # borehole
        has_borehole = haskey(meta.geometry, :borehole)
        if has_borehole
            # borehole_well = to_SSD_units(T, meta.geometry.borehole.gap_in_mm, u"mm")
            # borehole_height = crystal_height - borehole_well
            borehole_height = to_SSD_units(T, meta.geometry.borehole.depth_in_mm, u"mm")
            borehole_radius = to_SSD_units(T, meta.geometry.borehole.radius_in_mm, u"mm")
            semiconductor_geometry -= CSG.Cone{T}(CSG.ClosedPrimitive; 
                r = borehole_radius, 
                hZ = borehole_height / 2 + gap, 
                # origin = CartesianPoint{T}(0, 0, borehole_well + borehole_height/2 + gap)
                origin = CartesianPoint{T}(0, 0, crystal_height - borehole_height/2 + gap)
            )
        end

        # top outer taper
        # top_outer_taper_height = to_SSD_units(T, meta.geometry.taper.top.outer.height_in_mm, u"mm")
        top_outer_taper_height = to_SSD_units(T, meta.geometry.taper.top.height_in_mm, u"mm")
        # if :radius_in_mm in keys(meta.geometry.taper.top.outer)
            if :radius_in_mm in keys(meta.geometry.taper.top)
            # top_outer_taper_radius = to_SSD_units(T, meta.geometry.taper.top.outer.radius_in_mm, u"mm")
            top_outer_taper_radius = to_SSD_units(T, meta.geometry.taper.top.radius_in_mm, u"mm")
            top_outer_taper_angle = atan(top_outer_taper_radius, top_outer_taper_height)
        # elseif :angle_in_deg in keys(meta.geometry.taper.top.outer)
        elseif :angle_in_deg in keys(meta.geometry.taper.top)
            # top_outer_taper_angle = to_SSD_units(T, meta.geometry.taper.top.outer.angle_in_deg, u"°")
            top_outer_taper_angle = to_SSD_units(T, meta.geometry.taper.top.angle_in_deg, u"°")
            top_outer_taper_radius = top_outer_taper_height * tan(top_outer_taper_angle)
        else
            error("The top outer taper needs either radius_in_mm or angle_in_deg")
        end
        has_top_outer_taper = top_outer_taper_height > 0 && top_outer_taper_angle > 0
        if has_top_outer_taper
            r_center = crystal_radius - top_outer_taper_height * tan(top_outer_taper_angle) / 2
            hZ = top_outer_taper_height/2 + 1gap
            Δr = hZ * tan(top_outer_taper_angle)         
            r_in_bot = r_center + Δr
            r_in_top = r_center - Δr
            r_out = max(r_in_top, r_in_bot) + gap # ensure that r_out is always bigger as r_in
            r = ((r_in_bot, r_out),(r_in_top, r_out))
            semiconductor_geometry -= CSG.Cone{T}(CSG.ClosedPrimitive; 
                r = r,
                hZ = hZ, 
                origin = CartesianPoint{T}(0, 0, crystal_height - top_outer_taper_height/2)
            )
        end

        # top inner taper -> borehole taper
        # top_inner_taper_height = to_SSD_units(T, meta.geometry.taper.top.inner.height_in_mm, u"mm")
        top_inner_taper_height = to_SSD_units(T, meta.geometry.taper.borehole.height_in_mm, u"mm")
        # if :radius_in_mm in keys(meta.geometry.taper.top.inner)
            if :radius_in_mm in keys(meta.geometry.taper.borehole)
            # top_inner_taper_radius = to_SSD_units(T, meta.geometry.taper.top.inner.radius_in_mm, u"mm")
            top_inner_taper_radius = to_SSD_units(T, meta.geometry.taper.borehole.radius_in_mm, u"mm")
            top_inner_taper_angle = atan(top_inner_taper_radius, top_inner_taper_height)
        # elseif :angle_in_deg in keys(meta.geometry.taper.top.inner)
        elseif :angle_in_deg in keys(meta.geometry.taper.borehole)
            # top_inner_taper_angle = to_SSD_units(T, meta.geometry.taper.top.inner.angle_in_deg, u"°")
            top_inner_taper_angle = to_SSD_units(T, meta.geometry.taper.top.angle_in_deg, u"°")
            top_inner_taper_radius = top_inner_taper_height * tan(top_inner_taper_angle)
        else
            error("The top inner taper needs either radius_in_mm or angle_in_deg")
        end
        has_top_inner_taper = top_inner_taper_height > 0 && top_inner_taper_angle > 0
        if has_top_inner_taper && !has_borehole 
            error("A detector without a borehole cannot have a top inner taper.")
        end
        if has_top_inner_taper
            r_center = borehole_radius + top_inner_taper_height * tan(top_inner_taper_angle) / 2
            hZ = top_inner_taper_height/2
            Δr = hZ * tan(top_inner_taper_angle)         
            r_out_bot = r_center - Δr
	    r_out_top = r_center + Δr * (1 + 2*gap/hZ)
	    r = ((r_out_bot,), (r_out_top,))
            semiconductor_geometry -= CSG.Cone{T}(CSG.ClosedPrimitive; 
                r = r,
                hZ = hZ + gap, 
                origin = CartesianPoint{T}(0, 0, crystal_height - top_inner_taper_height/2 + gap)
            )
        end

        # bot outer taper
        # bot_outer_taper_height = to_SSD_units(T, meta.geometry.taper.bottom.outer.height_in_mm, u"mm")
        bot_outer_taper_height = to_SSD_units(T, meta.geometry.taper.bottom.height_in_mm, u"mm")
        # if :radius_in_mm in keys(meta.geometry.taper.bottom.outer)
            if :radius_in_mm in keys(meta.geometry.taper.bottom)
            # bot_outer_taper_radius = to_SSD_units(T, meta.geometry.taper.bottom.outer.radius_in_mm, u"mm")
            bot_outer_taper_radius = to_SSD_units(T, meta.geometry.taper.bottom.radius_in_mm, u"mm")
            bot_outer_taper_angle = atan(bot_outer_taper_radius, bot_outer_taper_height)
        # elseif :angle_in_deg in keys(meta.geometry.taper.bottom.outer)
        elseif :angle_in_deg in keys(meta.geometry.taper.bottom)
            # bot_outer_taper_angle = to_SSD_units(T, meta.geometry.taper.bottom.outer.angle_in_deg, u"°")
            bot_outer_taper_angle = to_SSD_units(T, meta.geometry.taper.bottom.angle_in_deg, u"°")
            bot_outer_taper_radius = bot_outer_taper_height * tan(bot_outer_taper_angle)
        else
            error("The bottom outer tape needs either radius_in_mm or angle_in_deg")
        end
        has_bot_outer_taper = bot_outer_taper_height > 0 && bot_outer_taper_angle > 0
        if has_bot_outer_taper
            r_center = crystal_radius - bot_outer_taper_height * tan(bot_outer_taper_angle) / 2
            hZ = bot_outer_taper_height/2 + 1gap
            Δr = hZ * tan(bot_outer_taper_angle)         
            r_in_bot = r_center - Δr
            r_in_top = r_center + Δr
            r_out = max(r_in_top, r_in_bot) + gap # ensure that r_out is always bigger as r_in
            r = ((r_in_bot, r_out),(r_in_top, r_out))
            semiconductor_geometry -= CSG.Cone{T}(CSG.ClosedPrimitive; 
                r = r,
                hZ = hZ, 
                origin = CartesianPoint{T}(0, 0, bot_outer_taper_height/2)
            )
        end

        # groove
        # groove_outer_radius = to_SSD_units(T, meta.geometry.groove.outer_radius_in_mm, u"mm")
        groove_outer_radius = to_SSD_units(T, meta.geometry.groove.radius_in_mm.outer, u"mm")
        groove_depth = to_SSD_units(T, meta.geometry.groove.depth_in_mm, u"mm")
        groove_inner_radius = to_SSD_units(T, meta.geometry.groove.radius_in_mm.inner, u"mm")
        # groove_width = to_SSD_units(T, meta.geometry.groove.width_in_mm, u"mm")
        # has_groove = groove_outer_radius > 0 && groove_depth > 0 && groove_width > 0
        has_groove = groove_outer_radius > 0 && groove_depth > 0 && groove_inner_radius > 0
        if has_groove
            hZ = groove_depth/2 + 1
            # r_in = groove_outer_radius - groove_width
            r_in = groove_inner_radius
            r_out = groove_outer_radius
            r = ((r_in, r_out), (r_in, r_out))
            semiconductor_geometry -= CSG.Cone{T}(CSG.ClosedPrimitive; 
                r = r, 
                hZ = groove_depth / 2 + gap, 
                origin = CartesianPoint{T}(0, 0, groove_depth/2 - gap)
            )
        end

        # bulletization -> Tommaso changed to model it as small taper, ignore for now
        # is_bulletized = !all(values(meta.geometry.bulletization) .== 0)
        # is_bulletized && @warn "Bulletization is not implemented yet, ignore for now."

        # extras
        haskey(meta.geometry, :extra) && @warn "Extras are not implemented yet, ignore for now."


        ### POINT CONTACT ###

        point_contact_geometry = begin
            # pc_radius = to_SSD_units(T, meta.geometry.contact.radius_in_mm, u"mm")
            pc_radius = to_SSD_units(T, meta.geometry.pp_contact.radius_in_mm, u"mm")
            # pc_depth = to_SSD_units(T, meta.geometry.contact.depth_in_mm, u"mm")
            pc_depth = to_SSD_units(T, meta.geometry.pp_contact.depth_in_mm, u"mm")

            CSG.Cone{T}(CSG.ClosedPrimitive; 
                r = pc_radius, 
                hZ = pc_depth / 2, 
                origin = CartesianPoint{T}(0, 0, pc_depth / 2)
            )
        end


        ### MANTLE CONTACT ###
        
        mantle_contact_geometry = begin # top plate
            top_plate = begin
                r = if !has_borehole 
                    !has_top_outer_taper ? crystal_radius : crystal_radius - top_outer_taper_radius
                else 
                    r_in = borehole_radius
                    r_out = crystal_radius
                    if has_top_inner_taper r_in += top_inner_taper_radius end
                    if has_top_outer_taper r_out -= top_outer_taper_radius end
                    ((r_in, r_out), (r_in, r_out))
                end
                CSG.Cone{T}(CSG.ClosedPrimitive; 
                    r = r, 
                    hZ = li_thickness / 2, 
                    origin = CartesianPoint{T}(0, 0, crystal_height - li_thickness / 2)
                )
            end
            mc_geometry = top_plate
            
            # borehole at outer taper
            if has_top_outer_taper
                Δr_li_thickness = li_thickness / cos(top_outer_taper_angle)
                hZ = top_outer_taper_height/2
                r_bot = crystal_radius 
                r_top = crystal_radius - top_outer_taper_radius
                r = ((r_bot - Δr_li_thickness, r_bot),(r_top - Δr_li_thickness, r_top))
                mc_geometry += CSG.Cone{T}(CSG.ClosedPrimitive; 
                    r = r,
                    hZ = hZ, 
                    origin = CartesianPoint{T}(0, 0, crystal_height - top_outer_taper_height/2)
                )
            end

            # contact in borehole
            if has_top_inner_taper
                Δr_li_thickness = li_thickness / cos(top_inner_taper_angle)
                hZ = top_inner_taper_height/2    
                r_bot = borehole_radius
                r_top = borehole_radius + top_inner_taper_radius
                r = ((r_bot, r_bot+Δr_li_thickness),(r_top, r_top+Δr_li_thickness))
                mc_geometry += CSG.Cone{T}(CSG.ClosedPrimitive; 
                    r = r,
                    hZ = hZ, 
                    origin = CartesianPoint{T}(0, 0, crystal_height - top_inner_taper_height/2)
                )

                hZ = (borehole_height - top_inner_taper_height) / 2
                r = ((borehole_radius, borehole_radius+Δr_li_thickness),(borehole_radius, borehole_radius+Δr_li_thickness))
                mc_geometry += CSG.Cone{T}(CSG.ClosedPrimitive; 
                    r = r,
                    hZ = hZ, 
                    origin = CartesianPoint{T}(0, 0, crystal_height - top_inner_taper_height - hZ)
                )
            elseif has_borehole # but no inner taper
                hZ = borehole_height / 2
                r = ((borehole_radius, borehole_radius+li_thickness),(borehole_radius, borehole_radius+li_thickness))
                mc_geometry += CSG.Cone{T}(CSG.ClosedPrimitive; 
                    r = r,
                    hZ = hZ, 
                    origin = CartesianPoint{T}(0, 0, crystal_height - hZ)
                )
            end

            if has_borehole
                r = borehole_radius + li_thickness
                mc_geometry += CSG.Cone{T}(CSG.ClosedPrimitive; 
                    r = r, 
                    hZ = li_thickness / 2, 
                    origin = CartesianPoint{T}(0, 0, crystal_height - borehole_height - li_thickness / 2)
                )
            end

            # outer surface of mantle contact
            begin
                r = ((crystal_radius-li_thickness, crystal_radius),(crystal_radius-li_thickness, crystal_radius))
                hZ = crystal_height
                if has_top_outer_taper hZ -= top_outer_taper_height end
                z_origin = hZ/2
                if has_bot_outer_taper 
                    hZ -= bot_outer_taper_height 
                    z_origin += bot_outer_taper_height/2
                end
                mc_geometry += CSG.Cone{T}(CSG.ClosedPrimitive; 
                    r = r, 
                    hZ = hZ / 2, 
                    origin = CartesianPoint{T}(0, 0, z_origin)
                )
            end

            # bottom outer taper contact
            if has_bot_outer_taper
                Δr_li_thickness = li_thickness / cos(bot_outer_taper_angle)
                hZ = bot_outer_taper_height/2
                r_bot = crystal_radius - bot_outer_taper_radius
                r_top = crystal_radius
                r = ((r_bot - Δr_li_thickness, r_bot),(r_top - Δr_li_thickness, r_top))
                mc_geometry += CSG.Cone{T}(CSG.ClosedPrimitive; 
                    r = r,
                    hZ = hZ, 
                    origin = CartesianPoint{T}(0, 0, hZ)
                )
            end  

            # bottom surface of mantle contact (only if it has a groove ?)
            if groove_outer_radius > 0
                r_in = groove_outer_radius 
                r_out = crystal_radius
                if has_bot_outer_taper r_out -= bot_outer_taper_radius end
                r = ((r_in, r_out), (r_in, r_out))
                mc_geometry += CSG.Cone{T}(CSG.ClosedPrimitive; 
                    r = r, 
                    hZ = li_thickness / 2, 
                    origin = CartesianPoint{T}(0, 0, li_thickness / 2)
                )
            end

            mc_geometry
        end

        temperature = T(env.crystal_temperature)
        material = SolidStateDetectors.material_properties[:HPGe]
        
        # Impurity Model: Information are stored in `meta.production.impcc`
        # For now: Constant impurity density: 
        #   n-type: positive impurity density
        #   p-type: negative impurity density
        # Assume p-type
        constant_impurity_density = ustrip(uconvert(u"m^-3", T(-9*1e9) * u"cm^-3"))
        impurity_density_model = SolidStateDetectors.CylindricalImpurityDensity{T}(
            (0, 0, constant_impurity_density), # offsets
            (0, 0, 0)                          # linear slopes
        )

        # Charge Drift Model: 
        # Use example ADL charge drift model from SSD (Crystal axis <100> is at φ = 0):
        adl_charge_drift_config_file = joinpath(dirname(dirname(pathof(SolidStateDetectors))), 
            "examples/example_config_files/ADLChargeDriftModel/drift_velocity_config.yaml")
        charge_drift_model = SolidStateDetectors.ADLChargeDriftModel{T}(adl_charge_drift_config_file);

        semiconductor = SolidStateDetectors.Semiconductor(temperature, material, impurity_density_model, charge_drift_model, semiconductor_geometry)

        # operation_voltage = :op_voltage in keys(env) ? T(env.op_voltage) : T(meta.characterization.l200_site.recommended_voltage_in_V)
        operation_voltage = env.operating_voltage > 0 ? T(env.operating_voltage) : T(meta.characterization.l200_site.recommended_voltage_in_V)
        # ToDo: what if metadata has null

        point_contact = SolidStateDetectors.Contact( zero(T), material, 1, "Point Contact", point_contact_geometry )
        mantle_contact = SolidStateDetectors.Contact( operation_voltage, material, 2, "Mantle Contact", mantle_contact_geometry )

        semiconductor, (point_contact, mantle_contact)

        passives = missing # possible holding structure around the detector
        virtual_drift_volumes = missing
        SolidStateDetector{T}( meta.name, semiconductor, [point_contact, mantle_contact], passives, virtual_drift_volumes )
    end
