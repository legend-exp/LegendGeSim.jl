using Plots

## ---------- detector geometry
using SolidStateDetectors
using LegendGeSim 

## SSD IVC
simulation1 = Simulation{Float32}(SSD_examples[:InvertedCoax])
plot(simulation1.detector, size = (400, 400))

## public IVC metadata
ssd_conf = LegendGeSim.ssd_config("public_ivc_metadata.json", LegendGeSim.Environment(90,4000))
simulation6 = Simulation(SolidStateDetector{Float32}(ssd_conf))
plot6 = plot(simulation6.detector, size = (400, 400))


## ---------- simulation

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