# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function paczynski(u::AbstractFloat)
    absu = abs(u)
    if absu >= 1e5
        # You can verify that:
        #   u=1e5; (u*u + 2)/(u*sqrt(u*u + 4)) === 1.0
        return copysign(1.0, u)
    elseif absu <= 1e-8
        # You can verify that:
        #   u = 1e-8; (u*u + 2)/(u*sqrt(u*u + 4)) === inv(u)
        inv(u)
    else
        u2 = u*u
        # Using `hypot' in place of the square root would be an overkill since
        #   1e-8 < |u| < 1e5
        return (u2 + 2)/(u*sqrt(u2 + 4))
    end
end

"""
    paczynski(u) -> amplification

### Purpose ###

Calculate gravitational microlensing amplification of a point-like source by a
single point-like lens.

### Explanation ###

Return the [gravitational
microlensing](https://en.wikipedia.org/wiki/Gravitational_microlensing)
amplification of a point-like source by a single point-like lens, using
Paczyński formula

```math
A(u) = \\frac{u^2 + 2}{u\\sqrt{u^2 + 4}}
```

where \$u\$ is the projected distance between the lens and the source in units
of [Einstein radii](https://en.wikipedia.org/wiki/Einstein_radius).

In order to speed up calculations for extreme values of \$u\$, the following
asyntotic expressions for \$A(u)\$ are used:

```math
A(u) =
\\begin{cases}
 1/u & |u| \\ll 1 \\\\
 \\text{sgn}(u) & |u| \\gg 1
\\end{cases}
```

### Arguments ###

* `u`: projected distance between the lens and the source, in units of Einstein
  radii

### Output ###

The microlensing amplification for the given distance.

### Example ###

Calculate the microlensing amplification for \$u = 10^{-10}\$, \$10^{-1}\$,
\$1\$, \$10\$, \$10^{10}\$:

```jldoctest
julia> using AstroLib

julia> paczynski.([1e-10, 1e-1, 1, 10, 1e10])
5-element Vector{Float64}:
  1.0e10
 10.037461005722337
  1.3416407864998738
  1.0001922892047386
  1.0
```

### Notes ###

The expression of \$A(u)\$ of microlensing amplification has been given by
Bohdan Paczyński in

* Paczynski, B. 1986, ApJ, 304, 1.
  DOI:[10.1086/164140](http://dx.doi.org/10.1086/164140),
  Bibcode:[1986ApJ...304....1P](http://adsabs.harvard.edu/abs/1986ApJ...304....1P)

The same expression was actually found by Albert Einstein half a century
earlier:

* Einstein, A. 1936, Science, 84, 506.
  DOI:[10.1126/science.84.2188.506](http://dx.doi.org/10.1126/science.84.2188.506),
  Bibcode:[1936Sci....84..506E](http://adsabs.harvard.edu/abs/1936Sci....84..506E)
"""
paczynski(u::Real) = paczynski(float(u))
