# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

# TODO: give sense to "-0" ("-0" is bitwise equal to "0").
ten(degrees::T, minutes::T, seconds::T) where {T<:AbstractFloat} =
    copysign(1, degrees)*(abs(degrees) + minutes/60 + seconds/3600)

ten(d::Real, m::Real=0, s::Real=0) =
    ten(promote(float(d), float(m), float(s))...)

# TODO: improve performance, if possible.  There are a couple of slow tests to
# make parsing of the string work.
function ten(coordinate::AbstractString)
    lowercase(coordinate) == "missing" && return missing
    coord_strip = replace(uppercase(coordinate), r"[NSEW]" => "")
    split_coord = parse.(Float32, split(coord_strip, r"\s|:"))
    split_coord[2] /=60.0
    if length(split_coord) == 3
        split_coord[3] /=3600.0
    end
    conv = sum(abs.(split_coord))
    # N + E are positive | S + W are negative
    if split_coord[1] < 0 || occursin(r"[SW]", uppercase(coordinate))
        # negative
        return conv * -1
    else
        # positive
        return conv
    end
end

ten(itr) = ten(itr...)


"""
    ten(deg[, min, sec]) -> decimal
    ten("deg:min:sec") -> decimal

### Purpose ###

Converts a sexagesimal number or string to decimal.

### Explanation ###

`ten` is the inverse of the `sixty` function.

### Arguments ###

`ten` takes as argument either three scalars (`deg`, `min`, `sec`) or a string.
The string should have the form `"deg:min:sec"` or `"deg min sec"`.  Also any
iterable like `(deg, min, sec)` or `[deg, min, sec]` is accepted as argument.

If minutes and seconds are not specified they default to zero.

### Formatting requirements
- Coordinates as a `String` separated by spaces (`"11 43 41"`) or colons (`"11:43:41"`)
- Must use negative sign (`"-11 43.52"`) or single-letter cardinal direction (`"11 43.52W"`)
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)
- Can mix colons and spaces (although it's bad practice

### Output ###

The decimal conversion of the sexagesimal numbers provided is returned.

### Method ###

The formula used for the conversion is
```math
\\mathrm{sign}(\\mathrm{deg})·\\left(|\\mathrm{deg}| + \\frac{\\mathrm{min}}{60} + \\frac{\\mathrm{sec}}{3600}\\right)
```

### Example ###

```jldoctest
julia> using AstroLib

julia> ten(-0.0, 19, 47)
-0.3297222222222222

julia> ten("+5:14:58")
5.249444444444444

julia> ten("-10 26")
-10.433333333333334

julia> ten((-10, 26))
-10.433333333333334
```

### Notes ###

These functions cannot deal with `-0` (negative integer zero) in numeric input.
If it is important to give sense to negative zero, you can either make sure to
pass a floating point negative zero `-0.0` (this is the best option), or use
negative minutes and seconds, or non-integer negative degrees and minutes.
"""
function ten end
