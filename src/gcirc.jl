using SkyCoords
"""
AstroCoordinate represents a point in the sky with flexible input units
"""
struct AstroCoordinate{T<:AbstractFloat}
    coord::SkyCoords.ICRSCoords

    function AstroCoordinate(ra::T, dec::T, units::Symbol=:radians) where T<:AbstractFloat
        if units == :radians
            new{T}(SkyCoords.ICRSCoords(ra, dec))
        elseif units == :hours_degrees
            new{T}(SkyCoords.ICRSCoords(ra * π/12.0, deg2rad(dec)))
        elseif units == :degrees
            new{T}(SkyCoords.ICRSCoords(deg2rad(ra), deg2rad(dec)))
        else
            throw(ArgumentError("units must be :radians, :hours_degrees, or :degrees"))
        end
    end
end

"""
    gcirc(units, ra1, dec1, ra2, dec2) -> angular_distance

### Purpose ###

Computes rigorous great circle arc distances using modern Julia astronomical coordinates.

### Explanation ###

Input positions can be specified in radians, sexagesimal right ascension and declination, 
or degrees. Uses SkyCoords.jl for precise astronomical calculations.

### Arguments ###

* `units`: integer, can be either 0, or 1, or 2. Describes units of inputs and output:
    * 0: everything (input right ascensions and declinations, and output distance) 
         is radians
    * 1: right ascensions are in decimal hours, declinations in decimal degrees,
         output distance in arc seconds
    * 2: right ascensions and declinations are in degrees, output distance in arc
         seconds
* `ra1`: right ascension or longitude of point 1
* `dec1`: declination or latitude of point 1
* `ra2`: right ascension or longitude of point 2
* `dec2`: declination or latitude of point 2

Both `ra1` and `dec1`, and `ra2` and `dec2` can be given as 2-tuples `(ra1, dec1)` 
and `(ra2, dec2)`.

### Output ###

Angular distance on the sky between points 1 and 2, as an `AbstractFloat`. See
`units` argument above for the units.

### Method ###

Uses SkyCoords.jl's implementation of spherical geometry calculations, which provides
high-precision results for astronomical coordinate transformations and distance
measurements.

### Examples ###

```jldoctest
julia> using AstroLib

julia> gcirc(0, 120, -43, 175, +22)
1.590442261600714

julia> # Using tuple notation
julia> gcirc(1, (12.0, -43.0), (15.0, 22.0))
5823.795873591674  # Result in arcseconds

julia> # Using 
"""
function gcirc(units::Integer, ra1::T, dec1::T, ra2::T, dec2::T)::T where {T<:AbstractFloat}
    units_sym = if units == 0
        :radians
    elseif units == 1
        :hours_degrees
    elseif units == 2
        :degrees
    else
        throw(DomainError(units, "units must be 0 (radians), 1 (hours, degrees) or 2 (degrees)"))
    end
    
    coord1 = AstroCoordinate(ra1, dec1, units_sym)
    coord2 = AstroCoordinate(ra2, dec2, units_sym)
    
    distance = SkyCoords.separation(coord1.coord, coord2.coord)
    
    return units == 0 ? distance : rad2sec(distance)
end

# Type promotion for mixed Real inputs
function gcirc(units::Integer, ra1::Real, dec1::Real, ra2::Real, dec2::Real)
    T = promote_type(float(typeof(ra1)), float(typeof(dec1)), 
                    float(typeof(ra2)), float(typeof(dec2)))
    return gcirc(units, T(ra1), T(dec1), T(ra2), T(dec2))
end

# Tuple input handlers
function gcirc(units::Integer, radec1::Tuple{Real,Real}, ra2::Real, dec2::Real)
    T = promote_type(float(typeof(radec1[1])), float(typeof(radec1[2])), 
                    float(typeof(ra2)), float(typeof(dec2)))
    return gcirc(units, T(radec1[1]), T(radec1[2]), T(ra2), T(dec2))
end

function gcirc(units::Integer, ra1::Real, dec1::Real, radec2::Tuple{Real,Real})
    T = promote_type(float(typeof(ra1)), float(typeof(dec1)), 
                    float(typeof(radec2[1])), float(typeof(radec2[2])))
    return gcirc(units, T(ra1), T(dec1), T(radec2[1]), T(radec2[2]))
end

function gcirc(units::Integer, radec1::Tuple{Real,Real}, radec2::Tuple{Real,Real})
    T = promote_type(float(typeof(radec1[1])), float(typeof(radec1[2])), 
                    float(typeof(radec2[1])), float(typeof(radec2[2])))
    return gcirc(units, T(radec1[1]), T(radec1[2]), T(radec2[1]), T(radec2[2]))
end

# Updated broadcasting implementation
function Base.Broadcast.broadcasted(::typeof(gcirc), units::Integer, 
                                  ra1::AbstractArray, dec1::AbstractArray, 
                                  ra2::AbstractArray, dec2::AbstractArray)
    T = promote_type(float(eltype(ra1)), float(eltype(dec1)), 
                    float(eltype(ra2)), float(eltype(dec2)))
    return T[gcirc(units, r1, d1, r2, d2) for (r1, d1, r2, d2) in zip(ra1, dec1, ra2, dec2)]
end

function Base.Broadcast.broadcasted(::typeof(gcirc), units::Integer, 
                                  ra1::AbstractArray, dec1::AbstractArray, 
                                  ra2::Number, dec2::Number)
    T = promote_type(float(eltype(ra1)), float(eltype(dec1)), 
                    float(typeof(ra2)), float(typeof(dec2)))
    return T[gcirc(units, r1, d1, T(ra2), T(dec2)) for (r1, d1) in zip(ra1, dec1)]
end

function Base.Broadcast.broadcasted(::typeof(gcirc), units::Integer, 
                                  ra1::Number, dec1::Number,
                                  ra2::AbstractArray, dec2::AbstractArray)
    T = promote_type(float(typeof(ra1)), float(typeof(dec1)), 
                    float(eltype(ra2)), float(eltype(dec2)))
    return T[gcirc(units, T(ra1), T(dec1), r2, d2) for (r2, d2) in zip(ra2, dec2)]
end

# Helper function for converting radians to arcseconds
rad2sec(x::T) where T<:AbstractFloat = x * T(206264.8062470963)