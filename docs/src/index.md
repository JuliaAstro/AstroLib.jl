# AstroLib.jl

[AstroLib](https://github.com/JuliaAstro/AstroLib.jl) is a package of small generic routines useful above all in astronomical and astrophysical context, written in [Julia](https://github.com/julialang/julia.jl).

Included are also translations of some [IDL Astronomy User’s Library](http://idlastro.gsfc.nasa.gov/homepage.html) procedures, which are released under terms of [BSD-2-Clause License](http://idlastro.gsfc.nasa.gov/idlfaq.html#A14). AstroLib’s functions are not drop-in replacement of those procedures, Julia standard data types are often used (e.g., `DateTime` type instead of generic string for dates) and the syntax may slightly differ.

An extensive error testing suite ensures old fixed bugs will not be brought back by future changes.

## Installation

AstroLib is available for Julia 1.0 and later versions, and can be installed with [Julia](https://github.com/julialang/julia.jl)'s built-in package manager. In a Julia session run the commands

```julia-repl
julia> import Pkg
julia> Pkg.update()
julia> Pkg.add("AstroLib")
```

Older versions are also available for Julia 0.4-0.6.

Note that, in order to work, a few functions require external files, which are automatically downloaded when building the package. Should these files be missing for some reason, you will be able to load the package but some functions may not work properly. You can manually build the package with

```julia-repl
julia> Pkg.build("AstroLib")
```

## Usage

After installing the package, you can start using AstroLib with

```julia
using AstroLib
```

Many functions in `AstroLib.jl` are compatible with [Measurements.jl](https://github.com/giordano/Measurements.jl) package, which allows you to define quantities with uncertainty and propagate the error when performing calculations according to [propagation of uncertainty rules](https://en.wikipedia.org/wiki/Propagation_of_uncertainty). For example:

```jldoctest
julia> using AstroLib, Measurements

julia> mag2flux(12.54 ± 0.03)
3.499e-14 ± 9.7e-16
```

## How Can I Help?

`AstroLib.jl` is developed on [GitHub](https://github.com/juliaastro/AstroLib.jl). You can contribute to the project in a number of ways: by translating more routines from IDL Astronomy User’s Library, or providing brand-new functions, or even improving existing ones (make them faster and more precise). Also bug reports are encouraged.

## License

The `AstroLib.jl` package is licensed under the MIT “Expat” License. The original author is Mosè Giordano.

## Notes

This project is a work-in-progress, only few procedures have been translated so far. In addition, function syntax may change from time to time. Check [TODO.md](https://github.com/JuliaAstro/AstroLib.jl/blob/master/TODO.md) out to see how you can help. Volunteers are welcome!

## Documentation

Every function provided has detailed documentation that can be [accessed](http://docs.julialang.org/en/stable/manual/documentation/#accessing-documentation) at Julia REPL with

```julia-repl
julia> ?FunctionName
```

or with

```julia-repl
julia> @doc FunctionName
```

## Related Projects

This is not the only effort to bundle astronomical functions written in Julia language. Other packages useful for more specific purposes are available at [JuliaAstro](https://juliaastro.github.io/).

Because of this, some of IDL AstroLib’s utilities are not provided in `AstroLib.jl` as they are already present in other Julia packages. Here is a list of such utilities:

-   `aper`, see [Photometry.jl](https://github.com/juliaastro/Photometry.jl) package
-   `asinh`, already present in Julia with the same name
-   `cirrange`, it is equivalent to `mod(x, 360)`.  To restrict a number to the
    range `[0, 2pi)` use `mod2pi(x)`
-   `cosmo_param`, see [Cosmology.jl](https://github.com/JuliaAstro/Cosmology.jl) package
-   `galage`, see [Cosmology.jl](https://github.com/JuliaAstro/Cosmology.jl) package
-   `glactc_pm`, see [SkyCoords.jl](https://github.com/kbarbary/SkyCoords.jl) package
-   `glactc`, see [SkyCoords.jl](https://github.com/kbarbary/SkyCoords.jl) package
-   `jplephinterp`, see [JPLEphemeris.jl](https://github.com/helgee/JPLEphemeris.jl) package
-   `jplephread`, see [JPLEphemeris.jl](https://github.com/helgee/JPLEphemeris.jl) package
-   `jplephtest`, see [JPLEphemeris.jl](https://github.com/helgee/JPLEphemeris.jl) package
-   `lumdist`, see [Cosmology.jl](https://github.com/JuliaAstro/Cosmology.jl) package
-   `readcol`, use [readdlm](http://docs.julialang.org/en/stable/stdlib/io-network/#Base.readdlm), part of Julia `Base.DataFmt` module. This is not a complete replacement for `readcol` but most of the time it does-the-right-thing even without using any option (it automatically identifies string and numerical columns) and you do not need to manually specify a variable for each column

In addition, there are similar projects for Python ([Python AstroLib](http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/pyasl.html)) and R ([Astronomy Users Library](http://rpackages.ianhowson.com/cran/astrolibR/)).
