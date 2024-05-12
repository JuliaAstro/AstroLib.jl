# This file is a part of AstroLib.jl. License is MIT "Expat".

# dpdt gives the time rate of change of the mean orbital quantities
# dpdt elements taken from https://ssd.jpl.nasa.gov/txt/p_elem_t1.txt
const dpd = @SMatrix [ 0.00000037  0.00001906 -0.00594749 -0.12534081  0.16047689 149472.67411175;
                       0.00000390 -0.00004107 -0.00078890 -0.27769418  0.00268329  58517.81538729;
                       0.00000562 -0.00004392 -0.01294668  0.0         0.32327364  35999.37244981;
                       0.00001847  0.00007882 -0.00813131 -0.29257343  0.44441088  19140.30268499;
                      -0.00011607 -0.00013253 -0.00183714  0.20469106  0.21252668   3034.74612775;
                      -0.00125060 -0.00050991  0.00193609 -0.28867794 -0.41897216   1222.49362201;
                      -0.00196176 -0.00004397 -0.00242939  0.04240589  0.40805281    428.48202785;
                       0.00026291  0.00005105  0.00035372 -0.00508664 -0.32241464    218.45945325;
                      -0.00031596  0.00005170  0.00004818 -0.01183482 -0.04062942    145.20780515]

const record = Dict(1=>"mercury", 2=>"venus", 3=>"earth", 4=>"mars", 5=>"jupiter",
                  6=>"saturn", 7=>"uranus", 8=>"neptune", 9=>"pluto")

function _helio(jd::T, num::Integer, radians::Bool) where {T<:AbstractFloat}

    if num<1 || num>9
        throw(DomainError("Input should be an integer in the range 1:9 denoting planet number"))
    end
    t = (jd - J2000) / JULIANCENTURY
    body = record[num]
    dpdt = dpd .* t
    a = planets[body].axis/AU + dpdt[num, 1]
    eccen = planets[body].ecc + dpdt[num, 2]
    n = deg2rad(0.9856076686 / (a * sqrt(a) ))
    inc = deg2rad(planets[body].inc + dpdt[num, 3])
    along = deg2rad(planets[body].asc_long + dpdt[num, 4])
    plong = deg2rad(planets[body].per_long + dpdt[num, 5])
    mlong = deg2rad(planets[body].mean_long + dpdt[num, 6])
    m = mlong - plong
    E = kepler_solver(m, eccen)
    nu = trueanom(E, eccen)
    hrad = a * (1 - eccen * cos(E))
    hlong = mod2pi(nu + plong)
    hlat = asin(sin(hlong - along) * sin(inc))
    if !radians
        return hrad, rad2deg(hlong), rad2deg(hlat)
    end
    return hrad, hlong, hlat
end

"""
    helio(jd, list[, radians=true]) -> hrad, hlong, hlat

### Purpose ###

Compute heliocentric coordinates for the planets.

### Explanation ###

The mean orbital elements for epoch J2000 are used. These are derived
from a 250 yr least squares fit of the DE 200 planetary ephemeris to a
Keplerian orbit where each element is allowed to vary linearly with
time. Useful mainly for dates between 1800 and 2050, this solution fits the
terrestrial planet orbits to ~25'' or better, but achieves only ~600''
for Saturn.

### Arguments ###

* `jd`: julian date, scalar or vector
* `num`: integer denoting planet number, scalar or vector
  1 = Mercury, 2 = Venus, ... 9 = Pluto
* `radians`(optional): if this keyword is set to
  `true`, than the longitude and latitude output are in radians rather than degrees.

### Output ###

* `hrad`: the heliocentric radii, in astronomical units.
* `hlong`: the heliocentric (ecliptic) longitudes, in degrees.
* `hlat`: the heliocentric latitudes in degrees.

### Example ###

(1) Find heliocentric position of Venus on August 23, 2000

```jldoctest
julia> using AstroLib

julia> helio(jdcnv(2000,08,23,0), 2)
(0.7213758288364316, 198.39093251916148, 2.887355631705488)
```

(2) Find the current heliocentric positions of all the planets

```jldoctest
julia> using AstroLib

julia> helio.([jdcnv(1900)], 1:9)
9-element Array{Tuple{Float64,Float64,Float64},1}:
 (0.4207394142180803, 202.60972662618906, 3.0503005607270532)
 (0.7274605731764012, 344.5381482401048, -3.3924346961624785)
 (0.9832446886519147, 101.54969268801035, 0.012669354526696368)
 (1.4212659241051142, 287.8531100442217, -1.5754626002228043)
 (5.386813769590955, 235.91306092135062, 0.9131692817310215)
 (10.054339927304339, 268.04069870870387, 1.0851704598594278)
 (18.984683376211326, 250.0555468087738, 0.05297087029604253)
 (29.87722677219009, 87.07244903504716, -1.245060583142733)
 (46.9647515992327, 75.94692594417324, -9.576681044165511)
```
### Notes ###

This program is based on the two-body model and thus neglects
interactions between the planets.

The coordinates are given for equinox 2000 and *not* the equinox
of the supplied date.

Code of this function is based on IDL Astronomy User's Library.
"""
helio(jd::Real, num::Integer, radians::Bool=false) =
    _helio(float(jd), num, radians)

function helio(jd::AbstractVector{P}, num::AbstractVector{<:Real},
               radians::Bool = false) where {P<:Real}
    if length(jd) != length(num)
        throw(ArgumentError("jd and num vectors should be of the same length"))
    end
    typejd = float(P)
    hrad_out = similar(jd,  typejd)
    hlong_out = similar(jd,  typejd)
    hlat_out = similar(jd,  typejd)
    for i in eachindex(jd)
        hrad_out[i], hlong_out[i], hlat_out[i] =
        helio(jd[i], num[i], radians)
    end
    return hrad_out, hlong_out, hlat_out
end
