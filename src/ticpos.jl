# This file is a part of AstroLib.jl. License is MIT "Expat".

function ticpos(deglen::T, pixlen::T, ticsize::T) where {T <: AbstractFloat}
    minpix = deglen * 60 / pixlen
    incr = minpix * ticsize

    if incr >= 30
        units = "Degrees"
    elseif incr > 0.5
        units = "Arc Minutes"
    else
        units = "Arc Seconds"
    end

    if incr >= 120
        incr = T(4)
    elseif incr >= 60
        incr = T(2)
    elseif incr >= 30
        incr = T(1)
    elseif incr > 15
        incr = T(30)
    elseif incr >= 10
        incr = T(15)
    elseif incr >= 5
        incr = T(10)
    elseif incr >= 2
        incr = T(5)
    elseif incr >= 1
        incr = T(2)
    elseif incr > 0.5
        incr = T(1)
    elseif incr >= 0.25
        incr = T(30)
    elseif incr >= 0.16
        incr = T(15)
    elseif incr >= 0.08
        incr = T(10)
    elseif incr >= 0.04
        incr = T(5)
    elseif incr >= 0.02
        incr = T(2)
    else
        incr = T(1)
    end

    if units == "Degrees"
        minpix = minpix / 60
    elseif units == "Arc Seconds"
        minpix = minpix * 60
    end

    ticsize = incr / minpix
    return ticsize, incr, units
end

"""
    ticpos(deglen, pixlen, ticsize) -> ticsize, incr, units

### Purpose ###

Specify distance between tic marks for astronomical coordinate overlays.

### Explanation ###

User inputs number an approximate distance between tic marks,
and the axis length in degrees. `ticpos` will return a distance
between tic marks such that the separation is a round multiple
in arc seconds, arc minutes, or degrees.

### Arguments ###

* `deglen`: length of axis in degrees, positive scalar
* `pixlen`: length of axis in plotting units (pixels), postive scalar
* `ticsize`: distance between tic marks (pixels).  This value will be
   adjusted by `ticpos` such that the distance corresponds to a round
   multiple in the astronomical coordinate.

### Output ###

The 3-tuple `(ticsize, incr, units)`:

* `ticsize`: distance between tic marks (pixels), positive scalar
* `incr`: incremental value for tic marks in round units given
   by the `units` parameter
* `units`: string giving units of ticsize, either 'Arc Seconds',
  'Arc Minutes', or 'Degrees'

### Example ###

Suppose a 512 × 512 image array corresponds to 0.2 × 0.2 degrees on the sky.
A tic mark is desired in round angular units, approximately every 75 pixels.
Then

```jldoctest
julia> using AstroLib

julia> ticpos(0.2, 512, 75)
(85.33333333333333, 2.0, "Arc Minutes")
```

i.e. a good tic mark spacing is every 2 arc minutes, corresponding
to 85.333 pixels.

### Notes ###

All the arguments taken as input are assumed to be positive in nature.

Code of this function is based on IDL Astronomy User's Library.
"""
ticpos(deglen::Real, pixlen::Real, ticsize::Real) =
    ticpos(promote(float(deglen), float(pixlen), float(ticsize))...)
