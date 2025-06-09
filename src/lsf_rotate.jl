# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _lsf_rotate{T<:AbstractFloat}(Δ::T, vsini::T, ɛ::T)
    n_2 = floor(ceil(2*vsini/Δ)/2)
    vel = collect(-n_2:n_2)*Δ
    e1 = 2*(1.0 - ɛ)
    e2 = pi*ɛ/2.0
    e3 = pi*(1.0 - ɛ/3.0)
    x = vel/vsini
    x1 = abs.(1.0 - x.^2)
    return vel, (e1*sqrt.(x1) .+ e2*x1)/e3
end

"""
    lsf_rotate(delta_v, v_sin_i[, epsilon = 0.3]) -> velocity_grid, lsf

### Purpose ###

Create a 1-d convolution kernel to broaden a spectrum from a rotating star.

### Explanation ###

Can be used to derive the broadening effect (LSF, line spread function) due to
rotation on a synthetic stellar spectrum.  Assumes constant limb darkening
across the disk.

### Arguments ###

* `delta_v`: numeric scalar giving the step increment (in km/s) in the output
  rotation kernel
* `v_sin_i`: the rotational velocity projected along the line of sight (km/s)
* `epsilon` (optional numeric argument): numeric scalar giving the
  limb-darkening coefficient, default = 0.6 which is typical for photospheric
  lines.  The specific intensity \$I\$ at any angle \$\\theta\$ from the
  specific intensity \$I_{\\text{cen}}\$ at the center of the disk is given by:

  `` I = I_{\\text{cen}}\\cdot(1 - \\varepsilon\\cdot(1 - \\cos(\\theta)))``

### Output ###

The 2-tuple (`velocity_grid`, `lsf`):

* `velocity_grid`: vector of velocity grid with the same number of elements as
  `lsf` (see below)
* `lsf`: the convolution kernel vector for the specified rotational velocity.
  The number of points in `lsf` will be always be odd (the kernel is symmetric)
  and equal to either `ceil(2*v_sin_i/delta_v)` or `ceil(2*v_sin_i/delta_v) + 1`,
  whichever number is odd.  Elements of `lsf` will always be of type
  `AbstractFloat`.  To actually compute the broadening, the spectrum should be
  convolved with the rotational `lsf`

### Example ###

Plot the line spread function for a star rotating at 90 km/s in velocity space
every 3 km/s.  Use [PyPlot.jl](https://github.com/stevengj/PyPlot.jl) for
plotting.

``` julia
using PyPlot
plot(lsf_rotate(3, 90)...)
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
lsf_rotate(Δ::Real, v::Real, ɛ::Real=0.6) =
    _lsf_rotate(promote(float(Δ), float(v), float(ɛ))...)
