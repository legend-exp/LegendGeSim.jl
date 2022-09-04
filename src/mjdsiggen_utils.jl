function SolidStateDetectors.ElectricPotential(sim::SigGenSetup, ::Type{T} = Float32) where {T}
    E_pot, W_pot, E_abs, E_r, E_z = LegendGeSim.MJDSigGen.read_fields(sim);
    r_axis = (0:(sim.rlen - 1)) * sim.xtal_grid / 1000 # SSD is in SI Units: m
    z_axis = (0:(sim.zlen - 1)) * sim.xtal_grid / 1000 # SSD is in SI Units: m
    ax1 = SolidStateDetectors.DiscreteAxis{T, :r0, :fixed, ClosedInterval{T}}(r_axis[1]..r_axis[end], r_axis)
    ax2 = SolidStateDetectors.DiscreteAxis{T, :reflecting, :reflecting, ClosedInterval{T}}(zero(T)..zero(T), zeros(T, 1))
    ax3 = SolidStateDetectors.DiscreteAxis{T, :reflecting, :reflecting, ClosedInterval{T}}(z_axis[1]..z_axis[end], z_axis)
    axs = (ax1, ax2, ax3)
    grid = Grid{T, 3, SolidStateDetectors.Cylindrical, typeof(axs)}( axs )
    data = Array{T, 3}(undef, length(ax1), length(ax2), length(ax3));
    data[:] = E_pot'[:];
    ElectricPotential(data, grid)
end

function SolidStateDetectors.WeightingPotential(sim::SigGenSetup, ::Type{T} = Float32) where {T}
    E_pot, W_pot, E_abs, E_r, E_z = LegendGeSim.MJDSigGen.read_fields(sim);
    r_axis = (0:(sim.rlen - 1)) * sim.xtal_grid / 1000 # SSD is in SI Units: m
    z_axis = (0:(sim.zlen - 1)) * sim.xtal_grid / 1000 # SSD is in SI Units: m
    ax1 = SolidStateDetectors.DiscreteAxis{T, :r0, :fixed, ClosedInterval{T}}(r_axis[1]..r_axis[end], r_axis)
    ax2 = SolidStateDetectors.DiscreteAxis{T, :reflecting, :reflecting, ClosedInterval{T}}(zero(T)..zero(T), zeros(T, 1))
    ax3 = SolidStateDetectors.DiscreteAxis{T, :reflecting, :reflecting, ClosedInterval{T}}(z_axis[1]..z_axis[end], z_axis)
    axs = (ax1, ax2, ax3)
    grid = Grid{T, 3, SolidStateDetectors.Cylindrical, typeof(axs)}( axs )
    data = Array{T, 3}(undef, length(ax1), length(ax2), length(ax3));
    data[:] = W_pot'[:];
    WeightingPotential(data, grid)
end