# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _eqpole(l::T, b::T, southpole::Bool) where {T<:AbstractFloat}
    sgn = southpole ? -1 : 1
    l = deg2rad(sgn*l)
    b = deg2rad(sgn*b)
    r = 18 * sqrt(2 * (1 - sin(b))) * 3.53553391
    return r .* sincos(l)
end

"""
    eqpole(l, b[; southpole = false]) -> x, y

### Purpose ###

Convert right ascension \$l\$ and declination \$b\$ to coordinate \$(x, y)\$
using an equal-area polar projection.

### Explanation ###

The output \$x\$ and \$y\$ coordinates are scaled to be in the range \$[-90, 90]\$
and to go from equator to pole to equator.  Output map points can be
centered on the north pole or south pole.

### Arguments ###

* `l`: longitude, scalar or vector, in degrees
* `b`: latitude, same number of elements as right ascension, in degrees
* `southpole` (optional boolean keyword): keyword to indicate that the plot is
  to be centered on the south pole instead of the north pole.  Default is
  `false`.

### Output ###

The 2-tuple \$(x, y)\$:

* ``x`` coordinate, same number of elements as right ascension, normalized to be
  in the range ``[-90, 90]``.
* ``y`` coordinate, same number of elements as declination, normalized to be
  in the range ``[-90, 90]``.

### Example ###

```jldoctest
julia> using AstroLib

julia> eqpole(100, 35, southpole=true)
(-111.18287262822456, -19.604540237028665)

julia> eqpole(80, 19)
(72.78853915267848, 12.83458333897169)
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
eqpole(l::Real, b::Real; southpole::Bool=false) =
    _eqpole(promote(float(l), float(b))..., southpole)

function eqpole(l::AbstractArray{L}, b::AbstractArray{<:Real};
                southpole::Bool=false) where {L<:Real}
    if length(l) != length(b)
        throw(ArgumentError("l and b arrays should be of the same length"))
    end
    typel = float(L)
    x = similar(l, typel)
    y = similar(l, typel)
    for i in eachindex(l)
        x[i], y[i] = eqpole(l[i], b[i], southpole=southpole)
    end
    return x, y
end
