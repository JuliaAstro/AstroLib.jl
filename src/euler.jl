# This file is a part of AstroLib.jl. License is MIT "Expat".

function _euler(ai::T, bi::T, select::Integer, FK4::Bool, radians::Bool) where {T<:AbstractFloat}

    if select>6 || select<1
        throw(DomainError(
            "Input for coordinate transformation should be an integer in the range 1:6"))
    end

    if FK4
        psi = (0.57595865315, 4.9261918136, 0.0, 0.0, 0.11129056012, 4.7005372834)
        stheta = (0.88781538514, -0.88781538514,0.39788119938, -0.39788119938,
                  0.86766174755, -0.86766174755)
        ctheta = (0.46019978478, 0.46019978478, 0.91743694670, 0.9174369467,
                  0.49715499774, 0.49715499774)
        phi = (4.9261918136,  0.57595865315, 0.0, 0.0, 4.7005372834, 0.11129056012)
    else
        psi = (0.574770433, 4.9368292465, 0.0, 0.0, 0.11142137093, 4.71279419371)
        stheta = (0.88998808748, -0.88998808748, 0.39777715593, -0.39777715593,
                  0.86766622025, -0.86766622025)
        ctheta = (0.45598377618, 0.45598377618, 0.91748206207, 0.91748206207,
                  0.49714719172, 0.49714719172)
        phi = (4.9368292465, 0.574770433, 0.0, 0.0, 4.71279419371, 0.11142137093)
    end

    if !radians
        ai = deg2rad(ai)
        bi = deg2rad(bi)
    end
    sa, ca = sincos(ai - phi[select])
    sb, cb = sincos(bi)
    x = (cb * sa, cb * ca, sb)
    bo = ctheta[select]*x[3] - stheta[select]*x[1]
    ao = mod2pi(atan(ctheta[select]*x[1] + stheta[select]*x[3], x[2]) + psi[select])
    bo = asin(bo)

    if radians
        return (ao, bo)
    end
    return (rad2deg(ao), rad2deg(bo))
end

"""
    euler(ai, bi, select[, FK4=true, radians=true])

### Purpose ###

Transform between Galactic, celestial, and ecliptic coordinates.

### Explanation ###

The function is used by the astro procedure.

### Arguments ###

* `ai`: input longitude, scalar or vector.
* `bi`: input latitude, scalar or vector.
* `select` : integer input specifying type of coordinate
  transformation.
  SELECT   From          To     | SELECT   From       To
     1   RA-Dec (2000) Galactic |   4    Ecliptic   RA-Dec
     2   Galactic      RA-DEC   |   5    Ecliptic   Galactic
     3   RA-Dec        Ecliptic |   6    Galactic   Ecliptic

* `FK4` (optional boolean keyword) : if this keyword is set to `true`,
  then input and output celestial and ecliptic coordinates should be
  given in equinox B1950. When `false`, by default, they should be given in
  equinox J2000.
* `radians` (optional boolean keyword) : if this keyword is set to
  `true`, all input and output angles are in radians rather than degrees.

### Output ###

a 2-tuple `(ao, bo)`:

* `ao`: output longitude in degrees.
* `bo`: output latitude in degrees.

### Example ###

Find the Galactic coordinates of Cyg X-1 (ra=299.590315, dec=35.201604)

```jldoctest
julia> using AstroLib

julia> euler(299.590315, 35.201604, 1)
(71.33498957116959, 3.0668335310640984)
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
euler(ai::Real, bi::Real, select::Integer; FK4::Bool=false, radians::Bool=false) =
    _euler(promote(float(ai), float(bi))..., select, FK4, radians)

euler(aibi::Tuple{Real, Real}, select::Integer; FK4::Bool=false, radians::Bool=false) =
    euler(aibi[1], aibi[2], select, FK4=FK4, radians=radians)

function euler(ai::AbstractVector{R}, bi::AbstractVector{<:Real}, select::Integer;
               FK4::Bool=false, radians::Bool=false) where {R<:Real}
    if length(ai) != length(bi)
        throw(ArgumentError("ai and bi arrays should be of the same length"))
    end
    typeai = float(R)
    ai_out  = similar(ai,  typeai)
    bi_out = similar(bi, typeai)
    for i in eachindex(ai)
        ai_out[i], bi_out[i] = euler(ai[i], bi[i], select,
                                        FK4=FK4, radians=radians)
    end
    return ai_out, bi_out
end
