# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

using AstroLib
using Test
using Dates
using SkyCoords
using Unitful
using UnitfulAstro
#using Documenter

@testset "AstroLib" begin
    include("utils-tests.jl")
    include("misc-tests.jl")

    #DocMeta.setdocmeta!(AstroLib, :DocTestSetup, :(using AstroLib), recursive=true)
    #doctest(AstroLib)
end

# Dummy calls to "show" for new data types, just to increase code coverage.
show(devnull, AstroLib.planets["mercury"])
show(devnull, AstroLib.observatories["ca"])
show(devnull, AstroLib.observatories["vbo"])
