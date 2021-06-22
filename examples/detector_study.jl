using LegendGeSim 
using Plots

## read or launch simulation

detector = LegendGeSim.simulate_detector("data/detector_V02160A.json")

##
# active volume in cm^3
avol = LegendGeSim.get_active_volume(detector.point_types)
# capacitance in pF
cap = LegendGeSim.calculate_capacitance(detector)

## electric potential

plot(
    plot(detector.electric_potential, φ = 20), # initial electric potential (boundary conditions)
    plot(detector.point_types), # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
    plot(detector.ρ), # charge density distribution
    plot(detector.ϵ), # dielectric distribution
    layout = (1, 4), size = (1600, 500)
)

## electric field

plot(detector.electric_field, φ = 0.0, size = (350, 500))
LegendGeSim.plot_electric_fieldlines!(detector, φ = 0.0)

## weighting potential

plot(
    plot(detector.weighting_potentials[1]),
    plot(detector.weighting_potentials[2]),
    size = (900, 700)
)