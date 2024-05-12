# This file is a part of AstroLib.jl. License is MIT "Expat".

function _co_aberration(jd::T, ra::T, dec::T, eps::T) where {T<:AbstractFloat}
    t = (jd - J2000) / JULIANCENTURY
    if isnan(eps)
        eps = true_obliquity(jd)
    end
    sunlong = sunpos(jd, radians=true)[3]
    e = @evalpoly t 0.016708634 -0.000042037 -0.0000001267
    pe = @evalpoly t 102.93735 1.71946 0.00046
    sd, cd = sincos(deg2rad(dec))
    ce = cos(eps)
    te = tan(eps)
    sp, cp = sincos(deg2rad(pe))
    ss, cs = sincos(sunlong)
    sa, ca = sincos(deg2rad(ra))
    t1 = (cs*ce*(te*cd - sa*sd) + ca*sd*ss)
    t2 = (cp*ce*(te*cd - sa*sd) + ca*sd*sp)
    d_ra = 20.49552*(e*(ca*cp*ce + sa*sp) - ca*cs*ce - sa*ss)/cd
    d_dec = 20.49552*(e*t2 - t1)
    return d_ra, d_dec
end

"""
    co_aberration(jd, ra, dec[, eps=NaN]) -> d_ra, d_dec

### Purpose ###

Calculate changes to right ascension and declination due to the effect
of annual aberration

### Explanation ###

With reference to Meeus, Chapter 23

### Arguments ###

* `jd`: julian date, scalar or vector
* `ra`: right ascension in degrees, scalar or vector
* `dec`: declination in degrees, scalar or vector
* `eps` (optional): true obliquity of the ecliptic (in radians). It will be
  calculated if no argument is specified.

### Output ###

The 2-tuple `(d_ra, d_dec)`:

* `d_ra`: correction to right ascension due to aberration, in arc seconds
* `d_dec`: correction to declination due to aberration, in arc seconds

### Example ###

Compute the change in RA and Dec of Theta Persei (RA = 2h46m,11.331s, Dec = 49d20',54.5'')
due to aberration on 2028 Nov 13.19 TD

```jldoctest
julia> using AstroLib

julia> jd = jdcnv(2028,11,13,4, 56)
2.4620887055555554e6

julia> co_aberration(jd,ten(2,46,11.331)*15,ten(49,20,54.54))
(30.04404628365077, 6.699400463119431)
```

d_ra = 30.04404628365103'' (≈ 2.003s)
d_dec = 6.699400463118504''

### Notes ###

Code of this function is based on IDL Astronomy User's Library.

The output d_ra is *not* multiplied by cos(dec), so that
apparent_ra = ra + d_ra/3600.

These formula are from Meeus, Chapters 23.  Accuracy is much better than 1
arcsecond. The maximum deviation due to annual aberration is 20.49'' and occurs when the
Earth's velocity is perpendicular to the direction of the star.

This function calls [`true_obliquity`](@ref) and [`sunpos`](@ref).
"""
co_aberration(jd::Real, ra::Real, dec::Real, eps::Real=NaN) =
    _co_aberration(promote(float(jd), float(ra), float(dec), float(eps))...)

function co_aberration(jd::AbstractVector{R}, ra::AbstractVector{R},
                       dec::AbstractVector{R}, eps::Real=NaN) where {R<:Real}
    if !(length(jd) == length(ra) == length(dec))
        throw(ArgumentError("jd, ra and dec vectors should be of the same length"))
    end
    typejd = float(R)
    ra_out  = similar(ra,  typejd)
    dec_out = similar(dec, typejd)
    for i in eachindex(jd)
        ra_out[i], dec_out[i] = co_aberration(jd[i], ra[i], dec[i], eps)
    end
    return ra_out, dec_out
end
