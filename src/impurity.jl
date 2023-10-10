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


struct RadfordImpurityDensity{T} <: SolidStateDetectors.AbstractImpurityDensity{T}
    # a + b*z + c*exp((z-L)/tau) -> needs at least 4 points
    a::T
    b::T 
    c::T 
    tau::T
    L::T
    det_z0::T
end

function SolidStateDetectors.get_impurity_density(
    idm::RadfordImpurityDensity, pt::SolidStateDetectors.AbstractCoordinatePoint{T}
    )::T where {T}
    cpt = CartesianPoint(pt)
    z = cpt[3]

    # the function parameters are in crystal axis coordinates i.e. z = 0 is seed end, z = L crystal length 
    # -> convert to detector coordiantes where z = 0 corresponds to p+ contact i.e. z -> det_z0 - z
    -(idm.a .+ idm.b * (idm.det_z0 .- z) .+ idm.c * exp.((idm.det_z0 .- z .- idm.L)/idm.tau)) 

end


"""
    impurity_density_model(crystal_metadata)

PropDict -> Vector{Float64}

Fit the impurity measurements in the crystal metadata
with David Radford's empirical function a + b*z + c*exp((z-L)/tau)
where z is crystal coordinate from seed end (z=0) to crystal length(z=L)

The fit is performed in crystal metadata coordinates
i.e. distance in mm and impurity values in 1e9 e/cm^3
which determines the units of the parameters a,b,c,tau returned in the vector of Float32.
"""
function impurity_density_model(crystal_metadata::PropDict)
    # get impurity measurements
    T = Float32
    # distance from seed end in mm (crystal axis)
    dist = T.(crystal_metadata.impurity_measurements.distance_from_seed_end_mm)
    # measurements corr to dist
    imp_values = T.(crystal_metadata.impurity_measurements.value_in_1e9e_cm3) #  .* (1e6 * 1e9) # SSD in e/m^3

    # old: quadratic fit
    # convert points from crystal axis to detector axis, convert to SSD unit
    # dist_det = (det_z0 .- dist) ./ 1e3 # SSD in m
    # a, b, c = poly_fit(dist_det, imp_values, 2)
    # QuadraticImpurityDensity{Float32}(a, b, c)

    # Radford fit -> in the future need to make settings to choose?
    # fit in crystal metadata coordinates
    # sometimes might not have measurements at seed end? but last measurement can be taken as crystal length
    L = dist[end] # in mm

    @. radford(z, p) = p[1] + p[2] * z + p[3] * exp.((z-L)/p[4])

    ## starting values for fit
    # f(z = 0) = a
    # ToDo: what to do in the case no measurement at 0? should not happen in principle
    a = imp_values[1]
    # 0.04 good starting value in mm and 1e9e/cm^3
    b = 0.04 # * 1e3 * (1e6 * 1e9)
    # f(z = L) = a + bL + c
    c = (imp_values[end] - a - b*L)
    # tau = 20 good starting value in mm
    tau = 20 # / 1e3
    p0 = [a, b, c, tau]

    ### fit in 1e9 e/cm^3 VS mm
    fit = LsqFit.curve_fit(radford, dist, imp_values, p0)

    fit.param
end


function impurity_density_model(meta::PropDict, crystal_metadata_path::AbstractString, simulator::PSSimulator)

    # get crystal metadata
    crystal_name = first(meta.name, length(meta.name)-1)
    crystal_path = joinpath(crystal_metadata_path, crystal_name*".json")

    if !isdir(crystal_metadata_path) 
        @error "Crystal metadata path $crystal_metadata_path does not exist"
    end
    if !ispath(crystal_path)
        @error "Crystal $(meta.production.crystal) does not exist in path $crystal_metadata_path"
    end

    crystal_dict = propdict(crystal_path)
    
    # get fit parameters in crystal metadata coordinates i.e. e/cm^3 VS mm
    fit_param = impurity_density_model(crystal_dict)

    # construct SSD impurity density or create siggen impurity input file
    impurity_density_model(meta, crystal_dict, fit_param, simulator)
end


"""
Convert fit parameters from crystal metadata units (1e9 e/cm^3 VS mm) into SSD units (e/m^3 VS m)
and return SolidStateDetector impurity density type RadfordImpurityDensity to be later used
when constructing SSD object.
"""
function impurity_density_model(meta::PropDict, crystal_metadata::PropDict, fit_param::Vector{Float64}, ::SSDSimulator)    
    # convert to SSD units 
    a,b,c,tau = fit_param
    a = a * 1e9 * 1e6 # e/m^3 
    b = b * 1e9 * 1e6 * 1e3 # [e/m^3] / m 
    c = c * 1e9 * 1e6 # e/m^3
    tau = tau / 1e3 # m

    # crystal length
    T = Float32
    dist = T.(crystal_metadata.impurity_measurements.distance_from_seed_end_mm)
    # sometimes might not have measurements at seed end? but last measurement can be taken as crystal length
    L = dist[end] / 1e3 # m

    # position of detector axis Z=0 (p+ contact) of this detector (slice) in the crystal from seed end
    det_z0 = crystal_metadata.slices[Symbol(meta.production.slice)].detector_offset_in_mm / 1e3 # m   
    # fit parameters 

    RadfordImpurityDensity(T(a),T(b),T(c),T(tau),T(L),T(det_z0))
end


## ---------- fieldgen


"""
Convert fit parameters from crystal metadata units (1e9 e/cm^3 VS mm) into Siggen units (1e10 e/cm^3 VS mm)
and create a .dat file (unformatted stream of Float32) to be later used when calling fieldgen().
"""
function impurity_density_model(meta::PropDict, crystal_metadata::PropDict, fit_param::Vector{Float64}, ::SiggenSimulator)   

    a,b,c,tau = fit_param 

    T = Float32
    dist = T.(crystal_metadata.impurity_measurements.distance_from_seed_end_mm)
    # sometimes might not have measurements at seed end? but last measurement can be taken as crystal length
    L = dist[end] #mm

    # points in crystal axis from 0 to L with step of 1 mm
    # always put 200 points to make all .dat files same length
	imp = Array{Float32}(undef, 200, 1) # impurities in 10^10 cm^-3

    # length is L+1 including both ends
    Npoints = Int(ceil(L))+1 # step of 1 mm
    zpoints = LinRange(0, L, Npoints) 
    for i in 1:Npoints
        # divide by 10 to convert from 1e9 e/cm^3 to 1e10 e/cm^3
        # invert the order of the values because with siggen z=0 is tail (dirtiest part)
        # i.e. z = L in our system (Mirion measurements)
        imp[Npoints-i+1] = -(a + b * zpoints[i] + c * exp((zpoints[i] - L)/tau)) / 10.
    end
    for i in Npoints+1 : 200
        imp[i] = 0
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


#### OUTDATED, providing linear impurity profile with quadratic correction inside of siggen config
# """
#     fieldgen_impurity(meta)

# PropDict -> Vector{Float64}

# Obtain impurity profile parameters for the siggen config, namely
#     impurity_z0 = net impurity concentration at Z=0, in 1e10 e/cm3
#     impurity_gradient = net impurity gradient, in 1e10 e/cm4
#     impurity_quadratic = net impurity difference from linear, at z=L/2, in 1e10 e/cm3
# """
# function impurity_density_model(meta::PropDict, ::SiggenSimulator)
#     # values scaled by 10 because metadata is in 1e9 e/cm3 while fieldgen is in 1e10 e/cm3
#     imp_values = meta.production.impcc.array[Symbol("value_in_1e9e/cm3")] ./ 10        
#     dist = meta.production.impcc.array.dist_from_contact_in_mm ./ 10 # in cm

#     # !! test version - temp V08682B.json with full crystal impurity listed from the point of view of Z = 0 of the detector
#     # !! in the future planning to store crystal impurity in a separate file, and link it in the detector metadata
#     # !! together with the location of Z=0 in the crystal

#     # impurity curve - currently quadratic, David R. uses more complex linear+exp empirical function
#     f_quad = CurveFit.curve_fit(Polynomial, dist, imp_values, 2)
#     det_h = meta.geometry.height_in_mm / 10 # in cm

#     # !! --------- .dat file with unformatted Float32 with 1mm step
#     # step = 0.1 # step of 1mm in cm
#     # crystal_L = 14 #cm #!! temp hard coded 140mm crystal length for V08682B
#     # npoints = Int(crystal_L/step)
#     # x = LinRange(0,crystal_L,npoints)

#     # impfile = "configs/$(meta.det_name)_impurity.dat"
#     # open(impfile, "w") do io
#     #     for xp in x
#     #         y = Float32(f_quad(xp))
#     #         # David R: .dat = stream of Float32 unformatted
#     #         print(io, "$y ")
#     #     end
#     # end
#     # @info "impurity profile saved to $impfile"

#     # --------- coefficients based on quadratic fit
#     # net impurity concentration at Z=0, in 1e10 e/cm3
#     # "I make sure that the values at the bottom and top of the detector
#     # are “correct”. So I read them off from my curve for the overall profile
#     imp_z0 = f_quad(0)
#     # net impurity gradient, in 1e10 e/cm4
#     imp_grad = (f_quad(det_h)- imp_z0) / det_h

#     ## quadratic correction
#     # net impurity difference from linear, at z=L/2, in 1e10 e/cm3
#     imp_quad = f_quad(det_h / 2) - (imp_z0 + imp_grad*det_h/2)

#     # minus because p-type impurity results in negative charge after holes run away
#     -round.([imp_z0, imp_grad, imp_quad], digits=3)
# end
