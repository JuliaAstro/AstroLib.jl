# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

"""
    cirrange(number[, max]) -> restricted_number

### Purpose ###

Force a number into a given range \$[0, \\text{max})\$.

### Argument ###

* `number`: the number to modify.  Can be a scalar or an array.
* `max` (optional numerical argument): specify the extremum of the range \$[0, \\text{max})\$
 into which the number should be restricted.  If omitted, defaults to `360.0`.

### Output ###

The converted number or array of numbers, as `AbstractFloat`.

### Example ###

Restrict an array of numbers in the range \$[0, 2\\pi)\$ as if they are angles
expressed in radians:

``` julia
cirrange([4pi, 10, -5.23], 2.0*pi)
# => 3-element Array{Float64,1}:
#     0.0
#     3.71681
#     1.05319
```

### Notes ###

This function does not support the `radians` keyword like IDL implementation.
Use `2.0*pi` as second argument to restrict a number to the same interval.
"""
cirrange(number::Real, max::Real=360.0) =
    mod(promote(float(number), float(max))...)
