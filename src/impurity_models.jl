"""
    struct MJDFieldGenImpurityModel{T} <: SolidStateDetectors.AbstractImpurityDensity{T}

Impurity density for SolidStateDetectors for comparison to MJDFieldGen.
The definition and parameterization is based on the source code of MJDFieldGen.

In SSD, SI units are used. 
"""
struct MJDFieldGenImpurityModel{T} <: SolidStateDetectors.AbstractImpurityDensity{T}
    xtal_radius::T # radius of the detector
    xtal_length_half::T # half length (in z) of the detector
    # float impurity_z0;          // net impurity concentration at Z=0, in 1e10 e/cm3
    impurity_z0::T
    # float impurity_gradient;    // net impurity gradient, in 1e10 e/cm4
    impurity_gradient::T
    # float impurity_quadratic;   // net impurity difference from linear, at z=L/2, in 1e10 e/cm3
    impurity_quadratic::T
    # float impurity_surface;     // surface impurity of passivation layer, in 1e10 e/cm2
    # impurity_surface::T # This is handled in SSD via fix space charge density and not impurity densities
    # float impurity_radial_add;  // additive radial impurity at outside radius, in 1e10 e/cm3
    impurity_radial_add::T
    # float impurity_radial_mult; // multiplicative radial impurity at outside radius (neutral=1.0)
    impurity_radial_mult::T
    # float impurity_rpower;      // power for radial impurity increase with radius
    impurity_rpower::T
end

function MJDFieldGenImpurityModel{T}(s::SigGenSetup) where {T}
    if s.impurity_surface != 0
        @warn """FieldGen impurity density has a surface charge.
        This is handled in SSD via static charge densities and not the impurity density.
        Make sure this is correctly passed to SSD.
        """
    end
    MJDFieldGenImpurityModel{T}(
        ustrip(uconvert(u"m", (s.xtal_radius) * u"mm")), 
        ustrip(uconvert(u"m", (s.xtal_length/2) * u"mm")), 
        ustrip(uconvert(u"m^(-3)", s.impurity_z0 * u"1e10*cm^(-3)")),        
        ustrip(uconvert(u"m^(-4)", s.impurity_gradient * u"1e10*cm^(-4)")),  
        ustrip(uconvert(u"m^(-3)", s.impurity_quadratic * u"1e10*cm^(-3)")),  
        ustrip(uconvert(u"m^(-3)", s.impurity_radial_add * u"1e10*cm^(-3)")), 
        s.impurity_radial_mult,
        s.impurity_rpower,            
    )
end

function SolidStateDetectors.get_impurity_density(
    cdm::MJDFieldGenImpurityModel{T},
    pt::SolidStateDetectors.AbstractCoordinatePoint{T}
) where {T}
    cyl = CylindricalPoint(pt)
    r, z = cyl[1], cyl[3] 
    # From FieldGen Source code `mjd_fieldgen.c`:
    ρ_z = begin # z component
        # imp_z[i] = e_over_E * grid*grid / 4.0 *
        # (setup->impurity_z0 +
        #  setup->impurity_gradient * z * 0.1 +
        #  setup->impurity_quadratic *
        #  (1.0 - SQ(z - setup->xtal_length/2.0) / SQ(setup->xtal_length/2.0)));
        # }
        # #define SQ(x) ((x)*(x))
        cdm.impurity_z0 + 
        cdm.impurity_gradient * z + 
        cdm.impurity_quadratic * (T(1) - (z - cdm.xtal_length_half)^2 / cdm.xtal_length_half^2);
    end
    ρ_r_m = begin # radial mult. component
        # imp_rm = 1.0 + (setup->impurity_radial_mult - 1.0f) *
        #     pow((double) r / setup->xtal_radius, setup->impurity_rpower);
        1 + (cdm.impurity_radial_mult - 1) * (r / cdm.xtal_radius)^cdm.impurity_rpower
    end
    ρ_r_a = begin # radial additive component
        # imp_ra = setup->impurity_radial_add * e_over_E * grid*grid / 4.0 *
        #     pow((double) r / setup->xtal_radius, setup->impurity_rpower);
        cdm.impurity_radial_add * (r / cdm.xtal_radius)^cdm.impurity_rpower
    end
    # imp_z[i] * imp_rm + imp_ra;
    return ρ_z * ρ_r_m + ρ_r_a
end



function get_impurity_density_poly_from_metadata(::Type{T}, imp_dict::PropDict) where {T}
    zs_from_contact = T.(imp_dict.array.dist_from_contact_in_mm) # in mm
    imp_levels_at_z = T.(imp_dict.array[Symbol("value_in_1e9e/cm3")] ./ 10)  # in T(1e10) * u"cm^-3"
    return Polynomials.fit(zs_from_contact, imp_levels_at_z)
end

function determine_MJDFieldGenImpurityParameter_from_metadata(::Type{T}, detector_metadata::PropDict) where {T}
    imp_poly = get_impurity_density_poly_from_metadata(T, detector_metadata.production.impcc) # returns in 1e10/cm3
    R = T(detector_metadata.geometry.radius_in_mm) 
    L = T(detector_metadata.geometry.height_in_mm) 
    H = T(L/2)
    imp_level_at_0 = imp_poly(T(0))
    imp_level_at_H = imp_poly(H)
    imp_level_at_L = imp_poly(L)
    sub_poly = Polynomials.fit(T[0, H, L], T[imp_level_at_0, imp_level_at_H, imp_level_at_L])

    # the following parameters are determined by solving the equations of  
    # SolidStateDetectors.get_impurity_density(cdm::MJDFieldGenImpurityModel{T}, pt::SolidStateDetectors.AbstractCoordinatePoint{T})
    # for the z component
    impurity_z0 = imp_level_at_0
    impurity_quadratic = -sub_poly.coeffs[3] * H^2
    impurity_gradient = (sub_poly.coeffs[2] - 2*impurity_quadratic/H) * 10  # *10 since unit conversion inv(cm^3 * mm)

    return (;
        impurity_z0,
        impurity_gradient,
        impurity_quadratic,
        impurity_surface = T(0),
        impurity_radial_add = T(0),
        impurity_radial_mult = T(1),
        impurity_rpower = T(0)
    )
end

function MJDFieldGenImpurityModel{T}(
    mjd_imp_par::NamedTuple{(:impurity_z0, 
                             :impurity_gradient, 
                             :impurity_quadratic, 
                             :impurity_surface, 
                             :impurity_radial_add, 
                             :impurity_radial_mult, 
                             :impurity_rpower)},
    detector_length, # in mm
    detector_radius  # in mm
    ) where {T}
    if mjd_imp_par.impurity_surface != 0
        @warn """FieldGen impurity density has a surface charge.
        This is handled in SSD via static charge densities and not the impurity density.
        Make sure this is correctly passed to SSD.
        """
    end
    MJDFieldGenImpurityModel{T}(
        ustrip(uconvert(u"m", detector_radius * u"mm")), 
        ustrip(uconvert(u"m", (detector_length/2) * u"mm")), 
        ustrip(uconvert(u"m^(-3)", mjd_imp_par.impurity_z0 * u"1e10*cm^(-3)")),        
        ustrip(uconvert(u"m^(-4)", mjd_imp_par.impurity_gradient * u"1e10*cm^(-4)")),  
        ustrip(uconvert(u"m^(-3)", mjd_imp_par.impurity_quadratic * u"1e10*cm^(-3)")),  
        ustrip(uconvert(u"m^(-3)", mjd_imp_par.impurity_radial_add * u"1e10*cm^(-3)")), 
        mjd_imp_par.impurity_radial_mult,
        mjd_imp_par.impurity_rpower,            
    )
end
