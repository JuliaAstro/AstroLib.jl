# AstroLib

| **Documentation**                       | [**Package Evaluator**][pkgeval-link] | **Build Status**                          | **Code Coverage**               |
|:---------------------------------------:|:-------------------------------------:|:-----------------------------------------:|:-------------------------------:|
| [![][docs-stable-img]][docs-stable-url] | [![][pkg-0.5-img]][pkg-0.5-url]       | [![Build Status][travis-img]][travis-url] | [![][coveral-img]][coveral-url] |
| [![][docs-latest-img]][docs-latest-url] | [![][pkg-0.6-img]][pkg-0.6-url]       | [![Build Status][appvey-img]][appvey-url] | [![][codecov-img]][codecov-url] |

Goal
----

The aim of this project is to provide users with a set of small generic routines
useful above all in astronomical and astrophysical context, written in
[Julia](http://julialang.org/).

Included are also translations of some
[IDL Astronomy User's Library](http://idlastro.gsfc.nasa.gov/homepage.html)
procedures, which are released under terms of
[BSD-2-Clause License](http://idlastro.gsfc.nasa.gov/idlfaq.html#A14).
AstroLib.jl's functions are not drop-in replacement of those procedures, Julia
standard data types are often used (e.g., `DateTime` type instead of generic
string for dates) and the syntax may slightly differ.  Refer to the
documentation of the functions for details.

**Note:** this project is a work-in-progress, only few procedures have been
translated so far.  In addition, function syntax may change from time to time.
Check [TODO.md](https://github.com/JuliaAstro/AstroLib.jl/blob/master/TODO.md) out
to see how you can help.  Volunteers are welcome!

Installation
------------

`AstroLib.jl` is available for Julia 0.6 and later versions, and can be
installed with
[Julia built-in package manager](http://docs.julialang.org/en/stable/manual/packages/).
In a Julia session run the command

```julia
julia> Pkg.update()
julia> Pkg.add("AstroLib")
```

Older versions are also available for Julia 0.4 and 0.5.

Note that, in order to work, a few functions require external files, which are
automatically downloaded when building the package.  Should these files be
missing for some reason, you will be able to load the package but some functions
may not work properly.  You can manually build the package with

```julia
julia> Pkg.build("AstroLib")
```

Usage
-----

After installing the package, you can start using it with

```julia
using AstroLib
```

Many functions in `AstroLib.jl` are compatible with
[`Measurements.jl`](https://github.com/giordano/Measurements.jl) package, which
allows you to define quantities with uncertainty and propagate the error when
performing calculations according to
[propagation of uncertainty rules](https://en.wikipedia.org/wiki/Propagation_of_uncertainty).
For example:

``` julia
using AstroLib, Measurements
mag2flux(12.54 ± 0.03)
# => 3.499451670283562e-14 ± 9.669342299577655e-16
```

New Types
---------

### Observatory ###

`AstroLib.jl` defines a new `Observatory` type.  This can be used to define a
new object holding information about an observing site.  It is a
[composite type](http://docs.julialang.org/en/stable/manual/types/#composite-types)
whose fields are

* `name` (`AbstractString` type): the name of the site
* `latitude` (`Real` type): North-ward latitude of the site in degrees
* `longitude` (`Real` type): East-ward longitude of the site in degrees
* `altitude` (`Real` type): altitude of the site in meters
* `tz` (`Real` type): the number of hours of offset from UTC

The type constructor `Observatory` can be used to create a new `Observatory`
object.  Its syntax is

``` julia
Observatory(name, lat, long, alt, tz)
```

`name` should be a string; `lat`, `long`, and `tz` should be anything that can
be converted to a floating number with `ten` function; `alt` should be a real
number.

A predefined list of some observing sites is provided with
`AstroLib.observatories` constant.  It is a dictionary whose keys are the
abbreviated names of the observatories.  For example, you can access information
of the European Southern Observatory with

``` julia
julia> obs = AstroLib.observatories["eso"]
Observatory: European Southern Observatory
latitude:    -29.256666666666668°N
longitude:   -70.73°E
altitude:    2347.0 m
time zone:   UTC-4

julia> obs.longitude
-70.73
```

You can list all keys of the dictionary with

``` julia
keys(AstroLib.observatories)
```

Feel free to contribute new sites or adjust information of already present ones.

### Planet ###

The package provides `Planet` type to hold information about Solar System
planets.  Its fields are

* Designation:

	* `name`: the name

* Physical characteristics:

	* `radius`: mean radius in meters
	* `eqradius`: equatorial radius in meters
	* `polradius`: polar radius in meters
	* `mass`: mass in kilogram

* Orbital characteristics (epoch J2000):

	* `ecc`: eccentricity of the orbit
	* `axis`: semi-major axis of the orbit in meters
	* `period`: sidereal orbital period in seconds

The constructor has this syntax:

``` julia
Planet(name, radius, eqradius, polradius, mass, ecc, axis, period)
```

The list of Solar System planets, from Mercury to Pluto, is available with
`AstroLib.planets` dictionary.  The keys of this dictionary are the lowercase
names of the planets.  For example:

``` julia
julia> AstroLib.planets["mercury"]
Planet:            Mercury
mean radius:       2.4397e6 m
equatorial radius: 2.4397e6 m
polar radius:      2.4397e6 m
mass:              3.3011e23 kg
eccentricity:      0.20563069
semi-major axis:   5.790905e10 m
period:            5.790905e10 s

julia> AstroLib.planets["mars"].eqradius
3.3962e6

julia> AstroLib.planets["saturn"].mass
5.6834e25
```

Documentation
-------------

Every function provided has detailed documentation that can be
[accessed](http://docs.julialang.org/en/stable/manual/documentation/#accessing-documentation)
at Julia REPL with

``` julia
julia> ?FunctionName
```

or with

``` julia
julia> @doc FunctionName
```

Full documentation of all functions can be accessed at
https://astrolibjl.readthedocs.io/.  There you can find the complete list of
all functions provided by this library.  You can also download the manual in PDF
format from https://media.readthedocs.org/pdf/astrolibjl/latest/astrolibjl.pdf.

How Can I Help?
---------------

You can contribute to the project in number of ways.  You can translate more
routines from IDL Astronomy User's Library or provide brand-new functions and
improve existing ones.  Also bug reports are encouraged.

Related Projects
----------------

This is not the only effort to bundle astronomical functions written in Julia
language.  Other packages useful for more specific purposes are available at
https://juliaastro.github.io/.

In addition, there are similar projects for Python
([Python AstroLib](http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/pyasl.html))
and R
([Astronomy Users Library](http://rpackages.ianhowson.com/cran/astrolibR/)).

License
-------

The AstroLib.jl package is licensed under the MIT "Expat" License.  The original
author is Mosè Giordano.


[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://astrolibjl.readthedocs.io/en/latest/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://astrolibjl.readthedocs.io/en/stable/

[pkgeval-link]: http://pkg.julialang.org/?pkg=AstroLib

[pkg-0.5-img]: http://pkg.julialang.org/badges/AstroLib_0.5.svg
[pkg-0.5-url]: http://pkg.julialang.org/detail/AstroLib.html

[pkg-0.6-img]: http://pkg.julialang.org/badges/AstroLib_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/detail/AstroLib.html

[travis-img]: https://travis-ci.org/JuliaAstro/AstroLib.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JuliaAstro/AstroLib.jl

[appvey-img]: https://ci.appveyor.com/api/projects/status/jfa9e54lv92rqd3m?svg=true
[appvey-url]: https://ci.appveyor.com/project/giordano/astrolib-jl

[coveral-img]: https://coveralls.io/repos/github/JuliaAstro/AstroLib.jl/badge.svg?branch=master
[coveral-url]: https://coveralls.io/github/JuliaAstro/AstroLib.jl?branch=master

[codecov-img]: https://codecov.io/gh/JuliaAstro/AstroLib.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaAstro/AstroLib.jl
