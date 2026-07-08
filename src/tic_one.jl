# This file is a part of AstroLib.jl. License is MIT "Expat".

function tic_one(zmin::T, pixx::T, incr::T, ra::Bool) where {T <: AbstractFloat}
    mul = ra ? 4 : 60
    min1 = zmin * mul
    incra = abs(incr)
    rem = min1 % incra
    sign = min1 * incr

    if sign > 0
        tic1 = pixx - abs(rem) * pixx / incra
        min2 = (min1 + incr - rem) / mul
    else
        tic1 = abs(rem) * pixx / incra
        min2 = (min1 - rem) / mul
    end
    return min2, tic1
end

"""
    tic_one(zmin, pixx, incr[, ra=true]) -> min2, tic1

### Purpose ###

Determine the position of the first tic mark for astronomical images.

### Explanation ###

For use in labelling images with right ascension and declination axes.
This routine determines the position in pixels of the first tic.

### Arguments ###

* `zmin`: astronomical coordinate value at axis zero point (degrees
   or hours).
* `pixx`: distance in pixels between tic marks (usually obtained from [`tics`](@ref)).
* `incr` - increment in minutes for labels (usually an even number obtained
   from the procedure [`tics`](@ref)).
* `ra` (optional boolean keyword): if true, incremental value being entered
   is in minutes of time, else it is assumed that value is in else it's in minutes of arc.
   Default is false.

### Output ###

The 2 tuple `(min2, tic1)`:

* `min2`: astronomical coordinate value at first tic mark
* `tic1`: position in pixels of first tic mark

### Example ###

Suppose a declination axis has a value of 30.2345 degrees at its
zero point.  A tic mark is desired every 10 arc minutes, which
corresponds to 12.74 pixels, with increment for labels being 10 minutes.
Then

```jldoctest
julia> using AstroLib

julia> tic_one(30.2345, 12.74, 10)
(30.333333333333332, 7.554820000000081)
```

yields values of min2 ≈ 30.333 and tic1 ≈ 7.55482, i.e. the first tic
mark should be labeled 30 deg 20 minutes and be placed at pixel value
7.55482.

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
tic_one(zmin::Real, pixx::Real, incr::Real, ra::Bool = false) =
    tic_one(promote(float(zmin), float(pixx), float(incr))..., ra)
