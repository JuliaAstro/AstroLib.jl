# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function kepler_solver(_M::Real, e::Real)
    @assert 0 <= e <= 1 "eccentricity must be in the range [0, 1]"
    # M must be in the range [-pi, pi], see Markley (1995), page 2.
    M = rem2pi(_M, RoundNearest)
    T = float(promote_type(typeof(M), typeof(e)))
    if iszero(M) || iszero(e)
        return T(M)
    end
    pi2 = abs2(T(pi))
    # equation (20)
    α = (3 * pi2 + 8 * (pi2 - pi * abs(M)) / (5 * (1 + e)))/(pi2 - 6)
    # equation (5)
    d = 3 * (1 - e) + α * e
    # equation (9)
    q = 2 * α * d * (1 - e) - M * M
    # equation (10)
    r = 3 * α * d * (d - 1 + e) * M + M * M * M
    # equation (14)
    w = cbrt(abs2(abs(r) + sqrt(q * q * q + r * r)))
    # equation (15)
    E1 = (2 * r * w / @evalpoly(w, q * q, q, 1) + M)/d
    # equation (26) & equation (27)
    f2, f3 = e .* sincos(E1)
    # equation (21)
    f0 = E1 - f2 - M
    # equation (25)
    f1 = 1 - f3
    # equation (22)
    δ3 = -f0 / (f1 - f0 * f2 / (2 * f1))
    # equation (23)
    δ4 = -f0 / @evalpoly(δ3, f1, f2 / 2, f3 / 6)
    # equations (24) and (28)
    δ5 = -f0 / @evalpoly(δ4, f1, f2 / 2, f3 / 6, - f2 / 24)
    return E1 + δ5 # equation 29
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

```math
M(t) = E(t) - e \\sin E(t)
```

where \$M(t) = (t - t_{0})/P\$ is the mean anomaly, \$E(t)\$ the eccentric
anomaly, \$e\$ the eccentricity of the orbit, \$t_0\$ is the time of periapsis
passage, and \$P\$ is the period of the orbit.  Usually the eccentricity is
given and one wants to find the eccentric anomaly \$E(t)\$ at a specific time
\$t\$, so that also the mean anomaly \$M(t)\$ is known.

### Arguments ###

* `M`: mean anomaly.
* `e`: eccentricity, in the elliptic motion regime (\$0 \\leq e \\leq 1\$)

### Output ###

The eccentric anomaly \$E\$, restricted to the range \$[-π, π]\$.

### Method ###

Many different numerical methods exist to solve Kepler's equation.  This
function implements the algorithm proposed in Markley (1995) Celestial Mechanics
and Dynamical Astronomy, 63, 101
(DOI:[10.1007/BF00691917](http://dx.doi.org/10.1007/BF00691917)).  This method
is not iterative, requires only four transcendental function evaluations, and
has been proved to be fast and efficient over the entire range of elliptic
motion \$0 \\leq e \\leq 1\$.

### Example ###

- Find the eccentric anomaly for an orbit with eccentricity \$e = 0.7\$ and
  for \$M(t) = 8π/3\$.

  ```jldoctest
  julia> using AstroLib

  julia> ecc = 0.7;

  julia> E = kepler_solver(8pi/3, ecc)
  2.5085279492864223
  ```

- Plot the eccentric anomaly as a function of mean anomaly for eccentricity
  \$e = 0\$, \$0.5\$, \$0.9\$.  Recall that `kepler_solver` gives \$E ∈ [-π,
  π]\$, use `mod2pi` to have it in \$[0, 2\\pi]\$.  Use
  [PyPlot.jl](https://github.com/JuliaPlots/Plots.jl/) for plotting.

  ```julia
  using AstroLib, PyPlot
  M = range(0, stop=2pi, length=1001)[1:end-1];
  for ecc in (0, 0.5, 0.9)
      plot(M, mod2pi.(kepler_solver.(M, ecc)))
  end
  ```

### Notes ###

The true anomaly can be calculated with `trueanom` function.
"""
function kepler_solver end
