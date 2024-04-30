# This file is a part of AstroLib.jl. License is MIT "Expat".

function precess_cd(cd::AbstractMatrix{T}, epoch1::T, epoch2::T, crval_old::AbstractVector{T},
                    crval_new::AbstractVector{T}, FK4::Bool) where {T<:AbstractFloat}
    t = (epoch2 - epoch1) / 1000

    if FK4
        st = (epoch1 - 1900) / 1000
        c = sec2rad(t*(20046.85 - st*(85.33 + st*0.37) + t*(-42.67 - st*0.37 - t*41.8)))
    else
        st = (epoch1 - 2000) / 1000
        c = sec2rad(t*(20043.109 - st*(85.33 + st*0.217) + t*(-42.665 - st*0.217 - t*41.8)))
    end
    pole_ra = zero(T)
    pole_dec = T(90)

    if (epoch1 == 2000 && epoch2 == 1950) || (epoch1 == 1950 && epoch2 == 2000)
        pole_ra, pole_dec = bprecess(pole_ra, pole_dec)
    else
        pole_ra, pole_dec = precess(pole_ra, pole_dec, epoch1, epoch2, FK4=FK4)
    end
    sind1, cosd1 = sincos(deg2rad(crval_old[2]))
    sind2, cosd2 = sincos(deg2rad(crval_new[2]))
    sinra = sind(crval_new[1] - pole_ra)
    sin_c, cos_c = sincos(c)
    cosfi = (cos_c - sind1 * sind2) / (cosd1 * cosd2)
    sinfi = (abs(sin_c) * sinra) / cosd1
    return cd * SMatrix{2,2}(cosfi, -sinfi, sinfi, cosfi)
end

"""
    precess_cd(cd, epoch1, epoch2, crval_old, crval_new[, FK4=true]) -> cd

### Purpose ###

Precess the coordinate description matrix.

### Explanation ###

The coordinate matrix is precessed from epoch1 to epoch2.

### Arguments ###

* `cd`: 2×2 coordinate description matrix in degrees
* `epoch1`: original equinox of coordinates, scalar
* `epoch2`: equinox of precessed coordinates, scalar
* `crval_old`: 2 element vector containing right ascension and declination
  in degrees of the reference pixel in the original equinox
* `crval_new`: 2 element vector giving crval in the new equinox
* `FK4` (optional boolean keyword): if this keyword is set to `true`,
  then the precession constants are taken in the FK4 reference frame. When it
  is `false`, the default is the FK5 frame

### Output ###

* `cd`: coordinate description containing precessed values

### Example ###

```jldoctest
julia> using AstroLib

julia> precess_cd([20 60; 45 45], 1950, 2000, [34, 58], [12, 83])
2×2 Matrix{Float64}:
  48.8944  147.075
 110.188   110.365
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
This function should not be used for values more than 2.5 centuries from the year 1900.
This function calls [`sec2rad`](@ref), [`precess`](@ref) and [`bprecess`](@ref).
"""
precess_cd(cd::AbstractMatrix{<:Real}, epoch1::Real, epoch2::Real, crval_old::AbstractVector{<:Real},
           crval_new::AbstractVector{<:Real}, FK4::Bool=false) =
               precess_cd(float(cd), promote(float(epoch1), float(epoch2))...,
                          float(crval_old), float(crval_new), FK4)
