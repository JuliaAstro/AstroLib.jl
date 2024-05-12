# This file is a part of AstroLib.jl. License is MIT "Expat".

"""
    ordinal(num) -> result

### Purpose ###

Convert an integer to a correct English ordinal string.

### Explanation ###
The first four ordinal strings are "1st", "2nd", "3rd", "4th" ....

### Arguments ###

* `num`: number to be made ordinal. It should be of type int.

### Output ###

* `result`: ordinal string, such as '1st' '3rd '164th' '87th' etc

### Example ###

```jldoctest
julia> using AstroLib

julia> ordinal.(1:5)
5-element Array{String,1}:
 "1st"
 "2nd"
 "3rd"
 "4th"
 "5th"
```

### Notes ###

This function does not support float arguments, unlike the IDL implementation.
Code of this function is based on IDL Astronomy User's Library.
"""
function ordinal(num::Integer)
    a = num % 100
    if a == 11 || a == 12 || a == 13
        suffix = "th"
    else
        a = num % 10
        if a == 1
            suffix = "st"
        elseif a == 2
            suffix = "nd"
        elseif a == 3
            suffix = "rd"
        else
            suffix = "th"
        end
    end
    return string(num) * suffix
end
