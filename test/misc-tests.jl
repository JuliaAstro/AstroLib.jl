# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

@testset "ordinal" begin
    @test ordinal.([3, 32, 391, 2412, 1000000]) ==
        ["3rd", "32nd", "391st", "2412th", "1000000th"]
end

# Test rad2sec
@testset "rad2sec" begin
    @test @inferred(rad2sec(1)) ≈ 206264.80624709636
    @test @inferred(rad2sec(pi)) ≈ 648000.0
    @test @inferred(sec2rad(rad2sec(12.34))) ≈ 12.34
end

# Test sec2rad
@testset "sec2rad" begin
    @test @inferred(sec2rad(3600 * 30)) ≈ pi / 6
    @test @inferred(1 / sec2rad(1)) ≈ 206264.80624709636
    @test @inferred(rad2sec(sec2rad(56.78))) ≈ 56.78
end
