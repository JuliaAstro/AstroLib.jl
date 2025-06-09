# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

"""
    get_juldate() -> julian_days

### Purpose ###

Return the number of Julian days for current time.

### Explanation ###

Return for current time the number of Julian calendar days since epoch
`-4713-11-24T12:00:00` as a floating point.

### Example ###

``` julia
get_juldate()
daycnv(get_juldate())
```

### Notes ###
Use `jdcnv` to get the number of Julian days for a different date.
"""
get_juldate() = jdcnv(now(Dates.UTC))
