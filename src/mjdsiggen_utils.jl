function SolidStateDetectors.ElectricPotential(sim::SigGenSetup)
    E_pot, W_pot, E_abs, E_r, E_z = LegendGeSim.MJDSigGen.read_fields(sim);
    T = eltype(E_pot)
    r_axis = (0:(sim.rlen - 1)) * sim.xtal_grid
    z_axis = (0:(sim.zlen - 1)) * sim.xtal_grid
    ax1 = SolidStateDetectors.DiscreteAxis{T, :r0, :fixed, ClosedInterval{T}}(r_axis[1]..r_axis[end], r_axis)
    ax2 = SolidStateDetectors.DiscreteAxis{T, :reflecting, :reflecting, ClosedInterval{T}}(zero(T)..zero(T), zeros(T, 1))
    ax3 = SolidStateDetectors.DiscreteAxis{T, :reflecting, :reflecting, ClosedInterval{T}}(z_axis[1]..z_axis[end], z_axis)
    axs = (ax1, ax2, ax3)
    grid = Grid{T, 3, SolidStateDetectors.Cylindrical, typeof(axs)}( axs )
    data = Array{T, 3}(undef, length(ax1), length(ax2), length(ax3));
    data[:] = E_pot'[:];
    ElectricPotential(data, grid)
end

function SolidStateDetectors.WeightingPotential(sim::SigGenSetup)
    E_pot, W_pot, E_abs, E_r, E_z = LegendGeSim.MJDSigGen.read_fields(sim);
    T = eltype(W_pot)
    r_axis = (0:(sim.rlen - 1)) * sim.xtal_grid
    z_axis = (0:(sim.zlen - 1)) * sim.xtal_grid
    ax1 = SolidStateDetectors.DiscreteAxis{T, :r0, :fixed, ClosedInterval{T}}(r_axis[1]..r_axis[end], r_axis)
    ax2 = SolidStateDetectors.DiscreteAxis{T, :reflecting, :reflecting, ClosedInterval{T}}(zero(T)..zero(T), zeros(T, 1))
    ax3 = SolidStateDetectors.DiscreteAxis{T, :reflecting, :reflecting, ClosedInterval{T}}(z_axis[1]..z_axis[end], z_axis)
    axs = (ax1, ax2, ax3)
    grid = Grid{T, 3, SolidStateDetectors.Cylindrical, typeof(axs)}( axs )
    data = Array{T, 3}(undef, length(ax1), length(ax2), length(ax3));
    data[:] = W_pot'[:];
    WeightingPotential(data, grid)
end