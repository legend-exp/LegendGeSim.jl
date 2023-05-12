## ---------- SSD

struct QuadraticImpurityDensity{T} <: SolidStateDetectors.AbstractImpurityDensity{T}
    a::T
    b::T
    c::T
end

function SolidStateDetectors.get_impurity_density(
    idm::QuadraticImpurityDensity, pt::SolidStateDetectors.AbstractCoordinatePoint{T}
    )::T where {T}
    cpt = CartesianPoint(pt)
    z = cpt[3]

    -(idm.a + idm.b * z + idm.c * z^2)
end


# note: reading and fitting is now the same for SSD or siggen, only outcome different (save file for siggen)
function impurity_density_model(meta::PropDict, ::SSDSimulator)
    # get crystal corresponding to detector 
    # in the future change to accessing LegendMetadata
    crystal_no = meta.production.crystal # XXX
    crystal_path = 

    # get array of impurity measurements
    T = Float32
    imp_values = T.(meta.production.impcc.array[Symbol("value_in_1e9e/cm3")]).* (1e6 * 1e9) # SSD in e/m^3
    dist = T.(meta.production.impcc.array.dist_from_contact_in_mm) ./ 1e3 # SSD in m

    # fit
    a, b, c = poly_fit(dist, imp_values, 2)

    QuadraticImpurityDensity{Float32}(a, b, c)
end


## ---------- fieldgen

"""
    fieldgen_impurity(meta)

PropDict -> Vector{Float64}

Obtain impurity profile parameters for the siggen config, namely
    impurity_z0 = net impurity concentration at Z=0, in 1e10 e/cm3
    impurity_gradient = net impurity gradient, in 1e10 e/cm4
    impurity_quadratic = net impurity difference from linear, at z=L/2, in 1e10 e/cm3
"""
function impurity_density_model(meta::PropDict, ::SiggenSimulator)
    # values scaled by 10 because metadata is in 1e9 e/cm3 while fieldgen is in 1e10 e/cm3
    imp_values = meta.production.impcc.array[Symbol("value_in_1e9e/cm3")] ./ 10        
    dist = meta.production.impcc.array.dist_from_contact_in_mm ./ 10 # in cm

    # !! test version - temp V08682B.json with full crystal impurity listed from the point of view of Z = 0 of the detector
    # !! in the future planning to store crystal impurity in a separate file, and link it in the detector metadata
    # !! together with the location of Z=0 in the crystal

    # impurity curve - currently quadratic, David R. uses more complex linear+exp empirical function
    f_quad = CurveFit.curve_fit(Polynomial, dist, imp_values, 2)
    det_h = meta.geometry.height_in_mm / 10 # in cm

    # !! --------- .dat file with unformatted Float32 with 1mm step
    # step = 0.1 # step of 1mm in cm
    # crystal_L = 14 #cm #!! temp hard coded 140mm crystal length for V08682B
    # npoints = Int(crystal_L/step)
    # x = LinRange(0,crystal_L,npoints)

    # impfile = "configs/$(meta.det_name)_impurity.dat"
    # open(impfile, "w") do io
    #     for xp in x
    #         y = Float32(f_quad(xp))
    #         # David R: .dat = stream of Float32 unformatted
    #         print(io, "$y ")
    #     end
    # end
    # @info "impurity profile saved to $impfile"

    # --------- coefficients based on quadratic fit
    # net impurity concentration at Z=0, in 1e10 e/cm3
    # "I make sure that the values at the bottom and top of the detector
    # are “correct”. So I read them off from my curve for the overall profile
    imp_z0 = f_quad(0)
    # net impurity gradient, in 1e10 e/cm4
    imp_grad = (f_quad(det_h)- imp_z0) / det_h

    ## quadratic correction
    # net impurity difference from linear, at z=L/2, in 1e10 e/cm3
    imp_quad = f_quad(det_h / 2) - (imp_z0 + imp_grad*det_h/2)

    # minus because p-type impurity results in negative charge after holes run away
    -round.([imp_z0, imp_grad, imp_quad], digits=3)

    # !! HADES tests
    # V08682A
    # [-1.252, 0.039, 0.12]
    # V08682B
    # [-2.3, 0.238, 0.15]
    # V09372A
    # [-1.748, 0.08543, 0.35]
    # V09374A
    # [-1.482, 0.0509, 0.01]
    # V09724A
    # [-1.499, 0.0531, 0.015]
end
