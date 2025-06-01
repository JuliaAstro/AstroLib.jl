# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _kepler_solver{T<:AbstractFloat}(M::T, e::T)
    if 0 <= e <= 1
        # M is assumed to be in the range [-pi, pi], see Markley (1995), page 2.
        # First restrict it to [0, 2pi], then move values above pi to [-pi, 0].
        M = cirrange(M, 2.0*pi)
        M > pi && (M -= 2.0*pi)
        if M == 0
            return 0.0
        else
            # equation (20)
            α = (3.0*pi*pi + 1.6*pi*(pi - abs(M))/(1.0 + e))/(pi*pi - 6.0)
            # equation (5)
            d = 3.0*(1.0 - e) + α*e
            # equation (9)
            q = 2.0*α*d*(1.0 - e) - M*M
            # equation (10)
            r = 3.0*α*d*(d - 1.0 + e)*M + M*M*M
            # equation (14)
            w = (abs(r) + sqrt(q*q*q + r*r))^(2.0/3.0)
            # equation (15)
            E1 = (2.0*r*w/(w*w + w*q + q*q) + M)/d
            # equation (21)
            f0 = E1 - e*sin(E1) - M
            # equation (25)
            f1 = 1.0 - e*cos(E1)
            # equation (26)
            f2 = e*sin(E1)
            # equation (27)
            f3 = 1.0 - f1
            # equation (28)
            f4 = -f2
            # equation (22)
            δ3 = -f0/(f1 - f0*f2/(2.0*f1))
            # equation (23)
            δ4 = -f0/(f1 + δ3*f2/2.0 + δ3*δ3*f3/6.0)
            # equation (24)
            δ5 = -f0/(f1 + δ4*f2/2.0 + δ4*δ4*f3/6.0 + δ4*δ4*δ4*f4/24.0)
            return E1 + δ5 # equation 29
        end
    else
        error("eccentricity must be in the range [0, 1]")
    end
end

"""
    kepler_solver(M, e) -> E

### Purpose ###

Solve Kepler's equation in the elliptic motion regime (\$0 \\leq e \\leq 1\$)
and return eccentric anomaly \$E\$.

### Explanation ###

In order to find the position of a body in elliptic motion (e.g., in the
two-body problem) at a given time \$t\$, one has to solve the [Kepler's
equation](https://en.wikipedia.org/wiki/Kepler%27s_equation)

``M(t) = E(t) - e\\sin E(t)``

where \$M(t) = (t - t_{0})/P\$ is the mean anomaly, \$E(t)\$ the eccentric
anomaly, \$e\$ the eccentricity of the orbit, \$t_0\$ is the time of periapsis
passage, and \$P\$ is the period of the orbit.  Usually the eccentricity is
given and one wants to find the eccentric anomaly \$E(t)\$ at a specific time
\$t\$, so that also the mean anomaly \$M(t)\$ is known.

### Arguments ###

* `M`: mean anomaly.  This can be either a scalar or an array
* `e`: eccentricity, in the elliptic motion regime (\$0 \\leq e \\leq 1\$)

### Output ###

The eccentric anomaly \$E\$, restricted to the range \$[-\\pi, \\pi]\$.  If an
array of mean anomalies is provided in input, an array of the same length as `M`
is returned.

### Method ###

Many different numerical methods exist to solve Kepler's equation.  This
function implements the algorithm proposed in Markley (1995) Celestial Mechanics
and Dynamical Astronomy, 63, 101
(DOI:[10.1007/BF00691917](http://dx.doi.org/10.1007/BF00691917)).  This method
is not iterative, requires only four transcendental function evaluations, and
has been proved to be fast and efficient over the entire range of elliptic
motion \$0 \\leq e \\leq 1\$.

### Example ###

(1) Find the eccentric anomaly for an orbit with eccentricity \$e = 0.7\$ and
for \$M(t) = 8\\pi/3\$.

``` julia
ecc = 0.7;
E = kepler_solver(8pi/3, ecc)
# => 2.5085279492864223
```

(2) Plot the eccentric anomaly as a function of mean anomaly for eccentricity
\$e = 0\$, \$0.5\$, \$0.9\$.  Recall that `kepler_solver` gives \$E \\in [-\\pi,
\\pi]\$, use `cirrange` to have it in \$[0, 2\\pi]\$.  Use
[PyPlot.jl](https://github.com/stevengj/PyPlot.jl) for plotting.

``` julia
using PyPlot
M=linspace(0, 2pi, 1001)[1:end-1];
for ecc in (0, 0.5, 0.9); plot(M, cirrange(kepler_solver(M, ecc), 2pi)); end
```

### Notes ###

The true anomaly can be calculated with `trueanom` function.
"""
kepler_solver(M::Real, e::Real) =
    _kepler_solver(promote(float(M), float(e))...)
