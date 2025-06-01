# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _deredd{T<:AbstractFloat}(Eby::T, by::T, m1::T, c1::T, ub::T)
    Rm1 = -0.33
    Rc1 = 0.19
    Rub = 1.53
    Eby0 = Eby > 0.0 ? Eby : 0.0
    return by - Eby0, m1 - Rm1*Eby0, c1 - Rc1*Eby0, ub - Rub*Eby0
end

"""
    deredd(Eby, by, m1, c1, ub) -> by0, m0, c0, ub0

### Purpose ###

Deredden stellar Stromgren parameters given for a value of E(b-y)

### Arguments ###

* `Eby`: color index E(b-y), scalar (E(b-y) = 0.73*E(B-V))
* `by`: b-y color (observed)
* `m1`: Stromgren line blanketing parameter (observed)
* `c1`: Stromgren Balmer discontinuity parameter (observed)
* `ub`: u-b color (observed)

All arguments can be either scalars or arrays all of the same length.

### Output ###

The 4-tuple `(by0, m0, c0, ub0)`.

* `by0`: b-y color (dereddened)
* `m0`: line blanketing index (dereddened)
* `c0`: Balmer discontinuity parameter (dereddened)
* `ub0`: u-b color (dereddened)

These are scalars or arrays of the same length as the input arguments.

### Example ###

``` julia
deredd(0.5, 0.2, 1.0, 1.0, 0.1)
# => (-0.3,1.165,0.905,-0.665)
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
deredd(Eby::Real, by::Real, m1::Real, c1::Real, ub::Real) =
    _deredd(promote(float(Eby), float(by), float(m1), float(c1), float(ub))...)

function deredd{E<:Real, B<:Real, M<:Real, C<:Real, U<:Real}(Eby::AbstractArray{E},
                                                             by::AbstractArray{B},
                                                             m1::AbstractArray{M},
                                                             c1::AbstractArray{C},
                                                             ub::AbstractArray{U})
    @assert length(Eby) == length(by) == length(m1) == length(c1) == length(ub)
    typeeby = typeof(float(one(E)))
    by0 = similar(Eby, typeeby)
    m0  = similar(Eby, typeeby)
    c0  = similar(Eby, typeeby)
    ub0 = similar(Eby, typeeby)
    for i in eachindex(Eby)
        by0[i], m0[i], c0[i], ub0[i] = deredd(Eby[i], by[i], m1[i], c1[i], ub[i])
    end
    return by0, m0, c0, ub0
end
