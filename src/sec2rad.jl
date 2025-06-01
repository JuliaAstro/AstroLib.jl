# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

_sec2rad{T<:AbstractFloat}(sec::T) = deg2rad(sec/3600.0)

"""
    sec2rad(sec) -> radians

### Purpose ###

Convert from seconds to radians.

### Argument ###

* `sec`: number of seconds.  It can be either a scalar or a vector.

### Output ###

The number of radians corresponding to `sec`.  If `sec` is an array, an array of
the same length is returned.

### Example ###

``` julia
sec2rad(3600*30)
# => 0.5235987755982988
```

### Notes ###

Use `rad2sec` to convert radians to seconds.
"""
sec2rad(sec::Real) = _sec2rad(float(sec))
