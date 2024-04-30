# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function planck_wave(wavelength::T, temperature::T) where {T<:AbstractFloat}
    c1  = 3.741771790075259e-16 # = 2*pi*h*c*c
    c2  = 1.43877735382772e-2   # = h*c/k
    return c1/(wavelength^5*(expm1(c2/wavelength/temperature)))
end

"""
    planck_wave(wavelength, temperature) -> black_body_flux

### Purpose ###

Calculate the flux of a black body per unit wavelength.

### Explanation ###

Return the spectral radiance of a black body per unit wavelength using [Planck's
law](https://en.wikipedia.org/wiki/Planck%27s_law)

```math
B_λ(λ, T) = \\frac{2hc^2}{λ^5}\\frac{1}{e^{\\frac{hc}{λ k_\\mathrm{B}T}} - 1}
```

### Arguments ###

* `wavelength`: wavelength at which the flux is to be calculated, in meters.
* `temperature`: the equilibrium temperature of the black body, in Kelvin.

### Output ###

The spectral radiance of the black body, in units of W/(sr·m³).

### Example ###

Plot the spectrum of a black body in \$[0, 3]\$ µm at \$5000\$ K.  Use
[PyPlot.jl](https://github.com/JuliaPlots/Plots.jl/) for plotting.

```julia
using PyPlot
wavelength = range(0, stop=3e-6, length=1000);
temperature = ones(wavelength)*5000;
flux = planck_wave.(wavelength, temperature);
plot(wavelength, flux)
```

### Notes ###

`planck_freq` calculates the flux of a black body per unit frequency.

Code of this function is based on IDL Astronomy User's Library.
"""
planck_wave(w::Real, t::Real) = planck_wave(promote(float(w), float(t))...)
