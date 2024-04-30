# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

"""
    ydn2md(year, day) -> date

### Purpose ###

Convert from year and day number of year to a date.

### Explanation ###

Returns the date corresponding to the `day` of `year`.

### Arguments ###

* `year`: the year, as an integer.
* `day`: the day of `year`, as an integer.

### Output ###

The date, of `Date` type, of \$\\text{day} - 1\$ days after January 1st of
`year`.

### Example ###

Find the date of the 60th and 234th days of the year 2016.

```jldoctest
julia> using AstroLib

julia> ydn2md.(2016, [60, 234])
2-element Vector{Dates.Date}:
 2016-02-29
 2016-08-21
```

### Note ###

`ymd2dn` converts from a date to day of the year.
"""
ydn2md(year::Integer, day::Integer) =
    Dates.firstdayofyear(Date(year)) + Dates.Day(day - 1)
