# keV_unit = 1.0u"keV"; keV = typeof(keV_unit);
# MHz_unit = 1.0u"MHz"; MHz = typeof(MHz_unit);
# const freq_unit = u"MHz"
ns_unit = 1.0u"ns"; ns = typeof(ns_unit);
μs_unit = 1.0u"μs"; μs = typeof(μs_unit);

T = Float32
germanium_ionization_energy = T(2.95)u"eV"

"""
A shelf to put all the preamplifier parameters in
later to be filled from a configuration file
"""
@with_kw struct PreAmp
    "PreAmp exp decay time"
    τ_decay::μs = T(5)*u"μs"
    
    "PreAmp rise time"
    τ_rise::ns = T(15)*u"ns"

    # "sigma for electronic noise"
    # noise_σ::Real = uconvert(u"eV", T(3)u"keV") / germanium_ionization_energy
end


"""
    function PreAmp(elec_conf)

Construct a struct with PreAmp parameters based on given electronics config PropDict
"""
function construct_PreAmp(sim_conf::PropDict)
    preamp = PreAmp(
        τ_decay=T(sim_conf.preamp.t_decay)*u"μs",
        τ_rise=T(sim_conf.preamp.t_rise)*u"ns"
        # noise_σ=uconvert(u"eV", T(elec_conf.noise)u"keV") / germanium_ionization_energy
    )

    preamp
end