# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

# TODO: give sense to "-0" ("-0" is bitwise equal to "0").
ten(degrees::T, minutes::T, seconds::T) where {T<:AbstractFloat} =
    copysign(1, degrees)*(abs(degrees) + minutes/60 + seconds/3600)

ten(d::Real, m::Real=0, s::Real=0) =
    ten(promote(float(d), float(m), float(s))...)

# TODO: improve performance, if possible.  There are a couple of slow tests to
# make parsing of the string work.
function ten(strng::AbstractString)
    # Convert strings into numbers, empty strings into 0s.
    # Replace in the string multiple spaces or colons with a single space, strip leading
    # whitespaces, and split the resulting string using the space as separator.
    tmp = map(x-> x=="" ? 0.0 : parse(Float64, x),
              split(lstrip(replace(strng, r"(\:| )+" => s" ")), " "))
    # Concatenate "tmp" with 3 zeros, so that "angle" has at least 3 elements
    # also with an empty string in input.
    angle = vcat(tmp, zeros(Float64, 3))::Vector{Float64}
    ten(angle[1], angle[2], angle[3])
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

### Output ###

The decimal conversion of the sexagesimal numbers provided is returned.

### Method ###

The formula used for the conversion is

\$\$\\mathrm{sign}(\\mathrm{deg})·\\left(|\\mathrm{deg}| + \\frac{\\mathrm{min}}{60} + \\frac{\\mathrm{sec}}{3600}\\right)\$\$

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
