# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

"""
    juldate(date::DateTime) -> reduced_julia_days

### Purpose ###

Convert from calendar to Reduced Julian Days.

### Explanation ###

Julian Day Number is a count of days elapsed since Greenwich mean noon on 1
January 4713 B.C.  Julian Days are the number of Julian days followed by the
fraction of the day elapsed since the preceding noon.

This function takes the given `date` and returns the number of Julian calendar
days since epoch `1858-11-16T12:00:00` (Reduced Julian Days = Julian Days -
2400000).

### Argument ###

* `date`: date in Julian Calendar, UTC standard.  It can be either e single date
  or an array of dates.  Each element can be given in `DateTime` type or
  anything that can be converted to that type.  In the case of vectorial input,
  each element is considered as a date, so you cannot provide a date by parts.

### Output ###

The number of Reduced Julian Days is returned.  If `date` is an array, an array
of the same length as `date` is returned.

### Example ###

Get number of Reduced Julian Days at 2016-03-20T15:24:00.

``` julia
juldate(DateTime(2016, 03, 20, 15, 24))
# => 57468.14166666667
juldate(2016, 03, 20, 15, 24)
# => 57468.14166666667
juldate("2016-03-20T15:24")
# => 57468.14166666667
```

### Notes ###

Julian Calendar is assumed, thus before `1582-10-15T00:00:00` this function is
*not* the inverse of `daycnv`.  For the conversion proleptic Gregorian date to
number of Julian days, use `jdcnv`, which is the inverse of `daycnv`.

Code of this function is based on IDL Astronomy User's Library.
"""
# Before 1582-10-15 Dates.datetime2julian uses proleptic Gregorian Calendar,
# instead AstroLib's juldate uses Julian Calendar.  In addition, after that day
# using Dates.datetime2julian(date)-2.4e6 would give results not precise enough.
# For all these reasons, we use here the same algorithm as AstroLib's juldate.
function juldate(dt::DateTime)
    year, month, day = Dates.yearmonthday(dt)
    hours, minutes, seconds, milliseconds =
        Dates.hour(dt), Dates.minute(dt),
        Dates.second(dt), Dates.millisecond(dt)

    if year < 0
        year += 1
    else
        year == 0 && error("There is no year zero in Julian Calendar")
    end

    day += hours/24 + minutes/1440 + seconds/86400 + milliseconds/86400000

    # Do not include leap year in January and February.
    if month < 3
        year -= 1
        month += 12
    end

    jd = fld(year, 4) + 365*(year - 1860) + floor(306001//10000*(month + 1)) +
         day - 1055//10

    # Adjust for Gregorian Calendar, started on 1582-10-15 (= RJD -100830.5).
    if jd > -100830.5
        a = div(year, 100)
        jd += 2 - a + fld(a, 4)
    end
    return jd
end

juldate(dt...) = juldate(DateTime(dt...))
