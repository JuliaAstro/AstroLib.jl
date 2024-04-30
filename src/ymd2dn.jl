# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

"""
    ymd2dn(date) -> number_of_days

### Purpose ###

Convert from a date to day of the year.

### Explanation ###

Returns the day of the year for `date` with January 1st being day 1.

### Arguments ###

* `date`: the date with `Date` type.  Can be a single date or an array of dates.

### Output ###

The day of the year for the given `date`.  If `date` is an array, returns an
array of days.

### Example ###

Find the days of the year for March 5 in the years 2015 and 2016 (this is a leap
year).

```jldoctest
julia> using AstroLib, Dates

julia> ymd2dn.([Date(2015, 3, 5), Date(2016, 3, 5)])
2-element Vector{Int64}:
 64
 65
```

### Note ###

`ydn2md` converts from year and day number of year to a date.
"""
const ymd2dn = Dates.dayofyear
