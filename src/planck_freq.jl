# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _planck_freq{T<:AbstractFloat}(frequency::T, temperature::T)
    c1  = 1.47449944028424e-50 # = 2*h/c*c
    c2  = 4.79924466221135e-11 # = h/k
    return c1*frequency^3/expm1(c2*frequency/temperature)
end

"""
    planck_freq(frequency, temperature) -> black_body_flux

### Purpose ###

Calculate the flux of a black body per unit frequency.

### Explanation ###

Return the spectral radiance of a black body per unit frequency using [Planck's
law](https://en.wikipedia.org/wiki/Planck%27s_law)

`` B_\\nu(\\nu, T) = \\frac{2h\\nu ^3}{c^2} \\frac{1}{e^\\frac{h\\nu}{k_\\mathrm{B}T} - 1} ``

### Arguments ###

* `frequency`: frequency at which the flux is to be calculated, in Hertz.
* `temperature`: the equilibrium temperature of the black body, in Kelvin.

Both arguments can be either scalar or arrays of the same length.

### Output ###

The spectral radiance of the black body, in units of W/(sr·m²·Hz).

### Example ###

Plot the spectrum of a black body in \$[10^{12}, 10^{15.4}]\$ Hz at \$8000\$ K.
Use [PyPlot.jl](https://github.com/stevengj/PyPlot.jl) for plotting.

``` julia
using PyPlot
frequency=logspace(12, 15.4, 1000);
temperature=ones(frequency)*8000;
flux=planck_freq(frequency, temperature);
plot(frequency, flux)
```

### Notes ###

`planck_wave` calculates the flux of a black body per unit wavelength.
"""
planck_freq(f::Real, t::Real) = _planck_freq(promote(float(f), float(t))...)
