# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

# Test cirrange
@test cirrange(12345) ≈ 105.0
@test cirrange.([3*e, 10, -86.95, 6*pi], 2*pi) ≈
    [1.8716601781975495, 3.7168146928204138, 1.0145943005142044, 0.0]

# Test rad2sec
@test rad2sec.([1]) ≈ [206264.80624709636]
@test rad2sec(pi) ≈ 648000.0
@test sec2rad(rad2sec(12.34)) ≈ 12.34

# Test sec2rad
@test sec2rad.([3600*30]) ≈ [pi/6]
@test 1/sec2rad(1) ≈ 206264.80624709636
@test rad2sec(sec2rad(56.78)) ≈ 56.78
