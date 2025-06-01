# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _flux2mag{T<:AbstractFloat}(flux::T, zero_point::T, ABwave::T)
    if isnan(ABwave)
        return -2.5*log10(flux) - zero_point
    else
        return -2.5*log10(flux) - 5.0*log10(float(ABwave)) - 2.406
    end
end

"""
    flux2mag(flux[, zero_point, ABwave=number]) -> magnitude

### Purpose ###

Convert from flux expressed in erg/(s cm² Å) to magnitudes.

### Explanation ###

This is the reverse of `mag2flux`.

### Arguments ###

* `flux`: the flux to be converted in magnitude, expressed in
  erg/(s cm² Å).  It can be either a scalar or an array.
* `zero_point`: scalar giving the zero point level of the magnitude.  If not
 supplied then defaults to 21.1 (Code et al 1976).  Ignored if the `ABwave`
 keyword is supplied
* `ABwave` (optional numeric keyword): wavelength scalar or vector in Angstroms.
 If supplied, then returns Oke AB magnitudes (Oke & Gunn 1983, ApJ, 266, 713;
 http://adsabs.harvard.edu/abs/1983ApJ...266..713O).

### Output ###

The magnitude.  It is of the same type, scalar or array, as `flux`.

If the `ABwave` keyword is set then magnitude is given by the expression

\$\$\\text{ABmag} = -2.5\\log_{10}(f) - 5\\log_{10}(\\text{ABwave}) - 2.406\$\$

Otherwise, magnitude is given by the expression

\$\$\\text{mag} = -2.5\\log_{10}(\\text{flux}) - \\text{zero point}\$\$

### Example ###

``` julia
flux2mag(5.2e-15)
# => 14.609991640913002
flux2mag(5.2e-15, 15)
# => 20.709991640913003
flux2mag(5.2e-15, ABwave=15)
# => 27.423535345634598
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
flux2mag(flux::Real, zero_point::Real=21.1; ABwave::Real=NaN) =
    _flux2mag(promote(float(flux), float(zero_point), float(ABwave))...)
