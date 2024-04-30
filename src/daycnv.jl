# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

"""
    daycnv(julian_days) -> DateTime

### Purpose ###

Converts Julian days number to Gregorian calendar dates.

### Explanation ###

Takes the number of Julian calendar days since epoch `-4713-11-24T12:00:00` and
returns the corresponding proleptic Gregorian Calendar date.

### Argument ###

* `julian_days`: Julian days number.

### Output ###

Proleptic Gregorian Calendar date, of type `DateTime`, corresponding to the
given Julian days number.

### Example ###

```jldoctest
julia> using AstroLib

julia> daycnv(2440000)
1968-05-23T12:00:00
```

### Notes ###

`jdcnv` is the inverse of this function.
"""
const daycnv = Dates.julian2datetime
