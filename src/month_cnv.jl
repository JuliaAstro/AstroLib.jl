# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

"""
    month_cnv(number[, shor=true, up=true, low=true]) -> month_name
    month_cnv(name) -> number

### Purpose ###

Convert between a month English name and  the equivalent number.

### Explanation ###

For example, converts from "January" to 1  or vice-versa.

### Arguments ###

The functions has two methods, one with numeric input (and three possible
boolean keywords) and the other one with string input.

Numeric input arguments:

* `number`: the number of the month to be converted to month name.
* `short` (optional boolean keyword): if true, the abbreviated (3-character)
  name of the month will be returned, e.g. "Apr" or "Oct".  Default is false.
* `up` (optional boolean keyword): if true, the name of the month will be all in
  upper case, e.g. "APRIL" or "OCTOBER".  Default is false.
* `low` (optional boolean keyword): if true, the name of the month will be all
  in lower case, e.g. "april" or "october".  Default is false.

String input argument:

* `name`: month name to be converted to month number.

### Output ###

The month name or month number, depending on the input.  For numeric input, the
format of the month name is influenced by the optional keywords.

### Example ###

```jldoctest
julia> using AstroLib

julia> month_cnv.(["janua", "SEP", "aUgUsT"])
3-element Vector{Int64}:
 1
 9
 8

julia> month_cnv.([2, 12, 6], short=true, low=true)
3-element Vector{String}:
 "feb"
 "dec"
 "jun"
```

"""
function month_cnv(number::Integer; short::Bool=false,
                   up::Bool=false, low::Bool=false)
    if short
        name = Dates.ENGLISH.months_abbr[number]
    else
        name = Dates.ENGLISH.months[number]
    end
    if up
        name = uppercase(name)
    elseif low
        name = lowercase(name)
    end
    return name
end

function month_cnv(name::AbstractString)
    name = strip(name)
    if length(name) >= 3
        output::Integer = getfield(Dates, Symbol(uppercasefirst(lowercase(name[1:3]))))
        return output
    else
        # Do the same as original AstroLib: return -1 for input string shorter
        # than 3.
        return -1
    end
end
