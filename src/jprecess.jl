# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

const Mjprec =
    SMatrix{6,6}(+0.9999256782,     +0.0111820610,     +0.0048579479,
                 -0.000551,         +0.238514,         -0.435623,
                 -0.0111820611,     +0.9999374784,     -0.0000271474,
                 -0.238565,         -0.002667,         +0.012254,
                 -0.0048579477,     -0.0000271765,     +0.9999881997,
                 +0.435739,         -0.008541,         +0.002117,
                 +0.00000242395018, +0.00000002710663, +0.00000001177656,
                 +0.99994704,       +0.01118251,       +0.00485767,
                 -0.00000002710663, +0.00000242397878, -0.00000000006582,
                 -0.01118251,       +0.99995883,       -0.00002714,
                 -0.00000001177656, -0.00000000006587, 0.00000242410173,
                 -0.00485767,       -0.00002718,       1.00000956)

# Note: IDL version of `jprecess' changes in-place "muradec", "parallax" and
# "radvel".  We don't do anything like this, but calculations are below,
# commented, in case someone is interested.
function _jprecess(ra::T, dec::T, parallax::T, radvel::T, epoch::T,
                   muradec::Vector{T}) where {T<:AbstractFloat}
    if length(muradec) != 2
        throw(DomainError("muradec must have length 2"))
    end
    sinra, cosra  = sincos(deg2rad(ra))
    sindec, cosdec = sincos(deg2rad(dec))
    if isfinite(epoch) && epoch != 1950
        A = A_precess .+ deg2rad.(A_dot_precess .* (epoch .- 1950) ./ 360000)
    else
        A = A_precess
    end
    r0 = SVector(cosra*cosdec,  sinra*cosdec,  sindec)
    r0_dot = SVector(-muradec[1]*sinra*cosdec - muradec[2]*cosra*sindec,
                     muradec[1]*cosra*cosdec - muradec[2]*sinra*sindec,
                     muradec[2]*cosdec) + 21.095 * radvel * parallax * r0
    r1 = r0 .- A .+ dot(r0, A) .* r0
    r1_dot = r0_dot .- A_dot_precess .+ dot(r0, A_dot_precess) .* r0
    R_1 = vcat(r1, r1_dot)
    R  = Mjprec*R_1
    if isfinite(epoch)
        t  = ((epoch - 1950) - 50.00021) / 100
        rr1 = @view(R[1:3]) .+ deg2rad.(@view(R[4:6]) .* t ./ 3600)
        x = rr1[1]
        y = rr1[2]
        z = rr1[3]
    else
        x = R[1]
        y = R[2]
        z = R[3]
        # x_dot = R[3]
        # y_dot= R[4]
        # z_dot = R[5]
    end
    rmag = norm((x, y, z))
    r2 = rmag*rmag
    dec2000 = asin(z / rmag)
    ra2000  = atan(y, x)
    # if isnan(epoch)
    #     muradec[1] = (x*y_dot - y*x_dot)/(x*x + y*y)
    #     muradec[2] = (z_dot*(x*x + y*y) - z*(x*x_dot + y*y_dot))/(r2*norm((x, y)))
    # end
    # if parallax > 0
    #     radvel   = (x*x_dot + y*y_dot + z*z_dot)/(21.095*parallax*rmag)
    #     parallax = parallax / rmag
    # end
    if ra2000 < 0
        ra2000 += 2 * T(pi)
    end
    ra2000  = rad2deg(ra2000)
    dec2000 = rad2deg(dec2000)
    return ra2000, dec2000
end

# Main interface.
jprecess(ra::Real, dec::Real, muradec::Vector{R};
         parallax::Real=0.0, radvel::Real=0.0) where {R<:Real} =
             _jprecess(promote(float(ra), float(dec), float(parallax),
                               float(radvel), NaN)...,
                       float(muradec))

jprecess(ra::Real, dec::Real, epoch::Real=1950.0) =
    _jprecess(promote(float(ra), float(dec), 0.0, 0.0, float(epoch))...,
              zeros(typeof(float(ra)), 2))

# Tuple arguments.
jprecess(radec::Tuple{Real,Real}, muradec::Vector{<:Real};
         parallax::Real=0.0, radvel::Real=0.0) =
             jprecess(radec..., muradec, parallax=parallax, radvel=radvel)

jprecess(radec::Tuple{Real,Real}, epoch::Real=1950.0) =
    jprecess(radec..., epoch)

# Vectorial arguments.
function jprecess(ra::AbstractArray{R}, dec::AbstractArray{<:Real},
                  muradec::AbstractArray{<:Real};
                  parallax::AbstractArray{<:Real}=zeros(R, length(ra)),
                  radvel::AbstractArray{<:Real}=zeros(R, length(ra))) where {R<:Real}
    if !(length(ra) == length(dec) == size(muradec)[2] == length(parallax) == length(radvel))
        # TODO write more helpful error message
        throw(ArgumentError(
            "ra, dec, muradec[:,2], parallax, and radvel must have the same length"))
    end
    typer = float(R)
    ra2000  = similar(ra, typer)
    dec2000 = similar(dec, typer)
    for i in eachindex(ra)
        ra2000[i], dec2000[i] = jprecess(ra[i], dec[i], muradec[:,i],
                                         parallax=parallax[i], radvel=radvel[i])
    end
    return ra2000, dec2000
end

function jprecess(ra::AbstractArray{R}, dec::AbstractArray{<:Real},
                  epoch::Real=1950.0) where {R<:Real}
    if length(ra) != length(dec)
        throw(ArgumentError("ra and dec arrays should be of the same length"))
    end
    typer = float(R)
    ra2000  = similar(ra, typer)
    dec2000 = similar(dec, typer)
    for i in eachindex(ra)
        ra2000[i], dec2000[i] = jprecess(ra[i], dec[i], epoch[i])
    end
    return ra2000, dec2000
end

"""
    jprecess(ra, dec[, epoch]) -> ra2000, dec2000
    jprecess(ra, dec, muradec[, parallax=parallax, radvel=radvel]) -> ra2000, dec2000

### Purpose ###

Precess positions from B1950.0 (FK4) to J2000.0 (FK5).

### Explanation ###

Calculate the mean place of a star at J2000.0 on the FK5 system from the mean
place at B1950.0 on the FK4 system.

`jprecess` function has two methods, one for each of the following cases:

* the proper motion is known and non-zero
* the proper motion is unknown or known to be exactly zero (i.e. extragalactic
  radio sources).  Better precision can be achieved in this case by inputting
  the epoch of the original observations.

### Arguments ###

The function has 2 methods.  The common mandatory arguments are:

* `ra`: input B1950 right ascension, in degrees.
* `dec`: input B1950 declination, in degrees.

The two methods have a different third argument (see "Explanation" section for
more details).  It can be one of the following:

* `muradec`: 2-element vector containing the proper motion in seconds of arc
  per tropical *century* in right ascension and declination.
* `epoch`: scalar giving epoch of original observations.

If none of these two arguments is provided (so `jprecess` is fed only with right
ascension and declination), it is assumed that proper motion is exactly zero and
`epoch = 1950`.

If it is used the method involving `muradec` argument, the following keywords
are available:

* `parallax` (optional numerical keyword): stellar parallax, in seconds of arc.
* `radvel` (optional numerical keyword): radial velocity in km/s.

Right ascension and declination can be passed as the 2-tuple `(ra, dec)`.  You
can also pass `ra`, `dec`, `parallax`, and `radvel` as arrays, all of the same
length N.  In that case, `muradec` should be a matrix 2×N.

### Output ###

The 2-tuple of right ascension and declination in 2000, in degrees, of input
coordinates is returned.  If `ra` and `dec` (and other possible optional
arguments) are arrays, the 2-tuple of arrays `(ra2000, dec2000)` of the same
length as the input coordinates is returned.

### Method ###

The algorithm is taken from the Explanatory Supplement to the Astronomical
Almanac 1992, page 184.  See also Aoki et al (1983), A&A, 128, 263.  URL:
http://adsabs.harvard.edu/abs/1983A%26A...128..263A.

### Example ###

The SAO catalogue gives the B1950 position and proper motion for the
star HD 119288.   Find the J2000 position.

* RA(1950) = 13h 39m 44.526s
* Dec(1950) = 8d 38' 28.63''
* Mu(RA) = -.0259 s/yr
* Mu(Dec) = -.093 ''/yr

```jldoctest
julia> using AstroLib

julia> muradec = 100 * [-15*0.0259, -0.093]; # convert to century proper motion

julia> ra = ten(13, 39, 44.526)*15;

julia> decl = ten(8, 38, 28.63);

julia> adstring(jprecess(ra, decl, muradec), precision=2)
" 13 42 12.740  +08 23 17.69"
```

### Notes ###

"When transferring individual observations, as opposed to catalog mean place,
the safest method is to tranform the observations back to the epoch of the
observation, on the FK4 system (or in the system that was used to to produce the
observed mean place), convert to the FK5 system, and transform to the the epoch
and equinox of J2000.0" -- from the Explanatory Supplement (1992), p. 180

`bprecess` performs the precession to B1950 coordinates.

Code of this function is based on IDL Astronomy User's Library.
"""
jprecess
