# This file is a part of LegendGeSim.jl, licensed under the MIT License (MIT).

using LegendGeSim
using Test

using LinearAlgebra, Random
# using Some, Other, Packages


@testset "utils" begin
    @testset "dummy_function" begin
        @test LegendGeSim.dummy_function(-4) == 16
        @test LegendGeSim.dummy_function(4) == 16
    end
end
