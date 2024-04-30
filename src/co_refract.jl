# This file is a part of AstroLib.jl. License is MIT "Expat".

function _co_refract(old_alt::T, altitude::T, pressure::T, temperature::T,
                    epsilon::T, to_observe::Bool) where {T<:AbstractFloat}

    if isnan(temperature)
        if altitude > 11000
            temperature = T(211.5)
        else
            temperature = 283 - 0.0065 * altitude
        end
    end

    if isnan(pressure)
        pressure = 1010 * (1 - 6.5/288000 * altitude) ^ 5.255
    end

    if !to_observe
        aout = old_alt - co_refract_forward(old_alt, pressure, temperature)
    else
        cur = old_alt + co_refract_forward(old_alt, pressure, temperature)
        while true
            last = cur
            cur = old_alt + co_refract_forward(last, pressure, temperature)
            if abs(last - cur) * 3600 < epsilon
                break
            end
        end
        aout = cur
    end
    return aout
end
"""
    co_refract(old_alt[, altitude=0, pressure=NaN, temperature=NaN,
               epsilon=0.25, to_observe=false]) -> aout

### Purpose ###

Calculate correction to altitude due to atmospheric refraction.

### Explanation ##

Because the index of refraction of air is not precisely 1.0, the atmosphere
bends all incoming light, making a star or other celestial object appear at
a slightly different altitude (or elevation) than it really is.  It is
important to understand the following definitions:

* Observed Altitude: The altitude that a star is seen to be, with a telescope.
  This is where it appears in the sky. This is should be always greater than the
  apparent altitude.

* Apparent Altitude: The altitude that a star would be at, if ~there were no
  atmosphere~ (sometimes called the "true" altitude). This is usually calculated
  from an object's celestial coordinates. Apparent altitude should always be
  smaller than the observed altitude.

Thus, for example, the Sun's apparent altitude when you see it right on the
horizon is actually -34 arcminutes.

This program uses a couple of simple formulae to estimate the effect for most
optical and radio wavelengths. Typically, you know your observed altitude
(from an observation), and want the apparent altitude. To go the other way,
this program uses an iterative approach.

### Arguments ###

* `old_alt`: observed altitude in degrees. If `to_observe` is set to true,
  this should be apparent altitude
* `altitude` (optional): the height of the observing location, in meters. This is
  only used to determine an approximate temperature and pressure, if these are not
  specified separately. Default is 0 i.e. sea level
* `pressure` (optional): the pressure at the observing location, in millibars.
  Default is NaN
* `temperature` (optional): the temperature at the observing location, in Kelvins.
  Default is NaN
* `epsilon` (optional): the accuracy to obtain, in arcseconds. If `to_observe` is
  true, then it will be calculated. Default is 0.25 arcseconds
* `to_observe` (optional boolean keyword): if set to true, it is assumed that `old_alt`
  has apparent altitude as its input and the observed altitude will be found

### Output ###

* `aout`: apparent altitude, in degrees. Observed altitude is returned if `to_observe`
  is set to true

### Example ###

The lower limb of the Sun is observed to have altitude of 0d 30'. Calculate the the true
(i.e. apparent) altitude of the Sun's lower limb using mean  conditions of air pressure and
temperature.

```jldoctest
julia> using AstroLib

julia> co_refract(0.5)
0.02584736873098442
```

### Notes ###

If altitude is set but the temperature or pressure is not, the program will make an
intelligent guess for the temperature and pressure.

#### Wavelength Dependence ####

This correction is 0 at zenith, about 1 arcminute at 45 degrees, and 34 arcminutes at the
horizon for optical wavelengths. The correction is non-negligible at all wavelengths,
but is not very easily calculable.
These formulae assume a wavelength of 550 nm, and will be accurate to about 4 arcseconds
for all visible wavelengths, for elevations of 10 degrees and higher.
Amazingly, they are also accurate for radio frequencies less than ~ 100 GHz.

#### References ####

* Meeus, Astronomical Algorithms, Chapter 15.
* Explanatory Supplement to the Astronomical Almanac, 1992.
* Methods of Experimental Physics, Vol 12 Part B, Astrophysics, Radio Telescopes, Chapter 2.5,
  "Refraction Effects in the Neutral Atmosphere", by R.K. Crane.

Code of this function is based on IDL Astronomy User's Library.
"""
co_refract(old_alt::Real, altitude::Real=0, pressure::Real=NaN, temperature::Real=NaN,
           epsilon::Real=0.25; to_observe::Bool=false) =
               _co_refract(promote(float(old_alt), float(altitude), float(pressure),
                          float(temperature), float(epsilon))..., to_observe)

"""
    co_refract_forward(alt, pre, temp) -> ref

### Purpose ###

A function used by [`co_refract`](@ref) to find apparent (or observed) altitude

### Arguments ###

* `alt`: the observed (or apparent) altitude, in degrees
* `pre`: pressure, in millibars
* `temp`: temperature, in Kelvins

### Output ###

* `ref`: the atmospheric refraction, in minutes of arc

### Notes ###

The atmospheric refraction is calculated by Saemundsson's formula

Code of this function is based on IDL Astronomy User's Library.
"""
function co_refract_forward(alt::Real, pre::Real, temp::Real)

    if alt < 15
        ref = 3.569*@evalpoly(alt, 0.1594, 0.0196, 0.00002)/@evalpoly(alt, 1, 0.505, 0.0845)
    else
        ref = 0.0166667 / tand((alt + 7.31 / (alt + 4.4)))
    end
    ref *= pre * 283 / (1010 * temp)
    return ref
end
