# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

_trueanom{T<:AbstractFloat}(E::T, e::T) =
    if e < 0
        # We don't need to check e > 1, because a DomainError is raised anyway.
        error("eccentricity must be in the range [0, 1]")
    else
        2.0*atan(sqrt((1.0 + e)/(1.0 - e))*tan(E/2.0))
    end

"""
    trueanom(E, e) -> true anomaly

### Purpose ###

Calculate true anomaly for a particle in elliptic orbit with eccentric anomaly
\$E\$ and eccentricity \$e\$.

### Explanation ###

In the two-body problem, once that the [Kepler's
equation](https://en.wikipedia.org/wiki/Kepler%27s_equation) is solved and
\$E(t)\$ is determined, the polar coordinates \$(r(t), \\theta(t))\$ of the body
at time \$t\$ in the elliptic orbit are given by

`` \\theta(t) = 2\\arctan \\left(\\sqrt{\\frac{1 + e}{1 - e}} \\tan\\frac{E(t)}{2} \\right)``

`` r(t) = \\frac{a(1 - e^{2})}{1 + e\\cos(\\theta(t) - \\theta_{0})}``

in which \$a\$ is the semi-major axis of the orbit, and \$\\theta_0\$ the value
of angular coordinate at time \$t = t_{0}\$.

### Arguments ###

* `E`: eccentric anomaly.  This can be either a scalar or an array
* `e`: eccentricity, in the elliptic motion regime (\$0 \\leq e \\leq 1\$)

### Output ###

The true anomaly.  If an array of eccentric anomalies is provided in input, an
array of the same length as `E` is returned.

### Example ###

Plot the true anomaly as a function of mean anomaly for eccentricity \$e = 0\$,
\$0.5\$, \$0.9\$.  Use [PyPlot.jl](https://github.com/stevengj/PyPlot.jl) for
plotting.

``` julia
using PyPlot
M=linspace(0, 2pi, 1001)[1:end-1];
for ecc in (0, 0.5, 0.9)
    E = kepler_solver(M, ecc);
    plot(M, cirrange(trueanom(E, ecc), 2pi))
end
```

### Notes ###

The eccentric anomaly can be calculated with `kepler_solver` function.
"""
trueanom(E::Real, e::Real) = _trueanom(promote(float(E), float(e))...)
