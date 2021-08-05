# Reference

## Data types

### Observatory

`AstroLib.jl` defines a new `Observatory` type. This can be used to define a new object holding information about an observing site. It is a [composite type] whose fields are

-   `name` (`String` type): the name of the site
-   `latitude` (`Float64` type): North-ward latitude of the site in degrees
-   `longitude` (`Float64` type): East-ward longitude of the site in degrees
-   `altitude` (`Float64` type): altitude of the site in meters
-   `tz` (`Float64` type): the number of hours of offset from UTC

The type constructor `Observatory` can be used to create a new `Observatory` object. Its syntax is

```julia
Observatory(name, lat, long, alt, tz)
```

`name` should be a string; `lat`, `long`, and `tz` should be anything that can be converted to a floating number with `ten` function; `alt` should be a real number.

A predefined list of some observing sites is provided with `AstroLib.observatories` constant. It is a dictionary whose keys are the abbreviated names of the observatories. For example, you can access information of the European Southern Observatory with

```julia
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

```julia
keys(AstroLib.observatories)
```

Feel free to contribute new sites or adjust information of already present ones.

### Planet

The package provides `Planet` type to hold information about Solar System planets. Its fields are

-   Designation:
    -   `name`: the name
-   Physical characteristics:
    -   `radius`: mean radius in meters
    -   `eqradius`: equatorial radius in meters
    -   `polradius`: polar radius in meters
    -   `mass`: mass in kilogram
-   Orbital characteristics (epoch J2000):
    -   `ecc`: eccentricity of the orbit
    -   `axis`: semi-major axis of the orbit in meters
    -   `period`: sidereal orbital period in seconds

The constructor has this syntax:

```julia
Planet(name, radius, eqradius, polradius, mass, ecc, axis, period)
```

The list of Solar System planets, from Mercury to Pluto, is available with `AstroLib.planets` dictionary. The keys of this dictionary are the lowercase names of the planets. For example:

```julia
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
5.6834e26
```

## Functions organized by category

### Coordinates and positions

[`adstring()`](@ref),
[`aitoff()`](@ref),
[`altaz2hadec()`](@ref),
[`baryvel()`](@ref),
[`bprecess()`](@ref),
[`co_aberration()`](@ref),
[`co_nutate()`](@ref),
[`co_refract()`](@ref),
[`eci2geo()`](@ref),
[`eq2hor()`](@ref),
[`eqpole()`](@ref),
[`euler()`](@ref),
[`gcirc()`](@ref),
[`geo2eci()`](@ref),
[`geo2geodetic()`](@ref),
[`geo2mag()`](@ref),
[`geodetic2geo()`](@ref),
[`hadec2altaz()`](@ref),
[`helio_rv()`](@ref),
[`helio()`](@ref),
[`hor2eq()`](@ref),
[`jprecess()`](@ref),
[`mag2geo()`](@ref),
[`mean_obliquity()`](@ref),
[`planet_coords()`](@ref),
[`polrec()`](@ref),
[`posang()`](@ref),
[`precess()`](@ref),
[`precess_cd()`](@ref),
[`precess_xyz()`](@ref),
[`premat()`](@ref),
[`radec()`](@ref),
[`recpol()`](@ref),
[`true_obliquity()`](@ref),
[`zenpos()`](@ref)

### Time and date

[`ct2lst()`](@ref),
[`daycnv()`](@ref),
[`get_date()`](@ref),
[`get_juldate()`](@ref),
[`helio_jd()`](@ref),
[`jdcnv()`](@ref),
[`juldate()`](@ref),
[`month_cnv()`](@ref),
[`nutate()`](@ref),
[`ydn2md()`](@ref),
[`ymd2dn()`](@ref)

### Moon and sun

[`moonpos()`](@ref),
[`mphase()`](@ref),
[`sunpos()`](@ref),
[`xyz()`](@ref)

### Utilities

[`airtovac()`](@ref),
[`calz_unred()`](@ref),
[`deredd()`](@ref),
[`flux2mag()`](@ref),
[`gal_uvw()`](@ref),
[`imf()`](@ref),
[`ismeuv()`](@ref),
[`kepler_solver()`](@ref),
[`lsf_rotate()`](@ref),
[`mag2flux()`](@ref),
[`paczynski()`](@ref),
[`planck_freq()`](@ref),
[`planck_wave()`](@ref),
[`rad2sec()`](@ref),
[`rhotheta()`](@ref),
[`sec2rad()`](@ref),
[`sixty()`](@ref),
[`sphdist()`](@ref),
[`ten()`](@ref),
[`tic_one()`](@ref),
[`ticpos()`](@ref),
[`tics()`](@ref),
[`trueanom()`](@ref),
[`uvbybeta()`](@ref),
[`vactoair()`](@ref)

### Miscellaneous (non-astronomy) functions

[`ordinal()`](@ref)

## Types and functions organized alphabetically

```@autodocs
Modules = [AstroLib]
```
