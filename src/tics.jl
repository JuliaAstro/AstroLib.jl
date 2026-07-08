function tics(radec_min::T, radec_max::T, numx::T, ticsize::T, ra::Bool) where {T <: AbstractFloat}

    numtics = numx / ticsize
    mul = ra ? T(4) : T(60)
    mins = abs(radec_min - radec_max) * mul
    rapix = (numx - 1) / mins
    incr = mins / numtics

    if incr >= 120
        incr = T(480)
    elseif incr >= 60
        incr = T(120)
    elseif incr >= 30
        incr = T(60)
    elseif incr >= 15
        incr = T(30)
    elseif incr >= 10
        incr = T(15)
    elseif incr >= 5
        incr = T(10)
    elseif incr >= 2
        incr = T(5)
    elseif incr >= 1
        incr = T(2)
    elseif incr >= 1 // 2
        incr = T(1)
    elseif incr >= 1 // 4
        incr = T(1 // 2)
    elseif incr >= 10 // 60
        incr = T(1 // 4)
    elseif incr >= 5 // 60
        incr = T(10 // 60)
    elseif incr >= 2 // 60
        incr = T(5 // 60)
    elseif incr >= 1 // 60
        incr = T(2 // 60)
    elseif incr >= 5 // 600
        incr = T(1 // 60)
    elseif incr >= 2 // 600
        incr = T(5 // 600)
    elseif incr >= 1 // 600
        incr = T(2 // 600)
    elseif incr >= 5 // 6000
        incr = T(1 // 600)
    elseif incr >= 2 // 6000
        incr = T(5 // 6000)
    elseif incr >= 1 // 6000
        incr = T(2 // 6000)
    else
        incr >= 0
        incr = T(1 // 6000)
    end

    ticsize = rapix * incr

    if ra && incr == 480
        incr = T(240)
    end
    if radec_min > radec_max
        incr = -incr
    end
    return ticsize, incr
end

"""
    tics(radec_min, radec_max, numx, ticsize[, ra=true]) -> ticsize, incr

### Purpose ###

Compute a nice increment between tic marks for astronomical images.

### Explanation ###

For use in labelling a displayed image with right ascension or declination
axes.  An approximate distance between tic marks is input, and a new value
is computed such that the distance between tic marks is in simple increments
of the tic label values.

### Arguements ###

* `radec_min` : minimum axis value (degrees).
* `radec_min` : maximum axis value (degrees).
* `numx` : number of pixels in x direction.
* `ticsize` : distance between tic marks (pixels).
* `ra` (optional boolean keyword): if true, incremental value would
  be in minutes of time. Default is false.

### Output ###

A 2-tuple `(ticsize, incr)`:
* `ticsize` : distance between tic marks (pixels).
* `incr` : incremental value for tic labels.  The format is dependent
  on the optional keyword. If true (i.e for right ascension), it's in minutes
  of time, else it's in minutes of arc (for declination).

### Example ###

```jldoctest
julia> using AstroLib

julia> tics(55, 60, 100.0, 1/2)
(0.66, 2.0)

julia> tics(30, 60, 12, 2, true)
(2.75, 30.0)
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
tics(
    radec_min::Real, radec_max::Real, numx::Real,
    ticsize::Real, ra::Bool = false
) = tics(
    promote(float(radec_min), float(radec_max), float(numx), float(ticsize))...,
    ra
)
