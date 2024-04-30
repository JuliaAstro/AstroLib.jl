# This file is a part of AstroLib.jl. License is MIT "Expat".

function ismeuv(wave::T, hcol::T, he1col::T, he2col::T, fano::Bool) where {T<:AbstractFloat}
    r = wave/911.75
    if r > 1
        return zero(T)
    end
    # minexp = -708.39642  from alog((machar(double=double).xmin)) where double = 1b in IDL
    z = sqrt(r/(1 - r))
    denom = -expm1(-2*T(π)*z)
    tauH = (3.44e-16 * r^4 *exp(-4*z*atan(1/z)) * hcol) / denom
    r *= 4
    tauHe2 = zero(T)
    if r < 1
        z = sqrt(r/(1 - r))
        denom = -4*expm1(-2*T(π)*z)
        tauHe2 = ((3.44e-16)*(r^4)*exp(-4*z*atan(1/z))*he2col)/denom
    end
    q  = (2.81, 2.51, 2.45, 2.44)
    fano_gamma = (2.64061e-03, 6.20116e-04, 2.56061e-04, 1.320159e-04)
    esubi = (4.421529414644497, 4.679309217126802, 4.738680412951545, 4.764345016358231)
    tauHe1 = zero(T)
    if wave < 503.97
        x = log10(wave)
        if wave < 46
            y = @evalpoly(x, -2.465188e+01, 4.354679, -3.553024, 5.573040,
                          -5.872938, 3.720797, -1.226919, 1.576657e-01)
        else
            y = @evalpoly(x, -2.953607e+01, 7.083061, 8.678646e-01, -1.221932,
                          4.052997e-02, 1.317109e-01, -3.265795e-02,
                          2.500933e-03)
            if fano
                episilon = 911.2671/wave
                for i = 1:4
                    x = 2*(episilon - esubi[i])/fano_gamma[i]
                    y += log10(((x - q[i])^2)/(1 + x^2))
                end
            end
        end
        tauHe1 = exp10(y)*he1col
    end
    return tauH + tauHe1 + tauHe2
end

"""
    ismeuv(wave, hcol[, he1col=hcol*0.1, he2col=0, fano=false]) -> tau

### Purpose ###

Compute the continuum interstellar EUV optical depth

### Explanation ###

The EUV optical depth is computed from the photoionization of hydrogen and helium.

### Arguments ###

* `wave`: wavelength value (in Angstroms). Useful range is 40 - 912 A;
  at shorter wavelength metal opacity should be considered, at longer wavelengths
  there is no photoionization.
* `hcol`: interstellar hydrogen column density in cm-2.
* `he1col` (optional): neutral helium column density in cm-2.
  Default is 0.1*hcol (10% of hydrogen column)
* `he2col` (optional): ionized helium column density in cm-2
  Default is 0.
* `fano` (optional boolean keyword): If this keyword is true, then the 4 strongest
  auto-ionizing resonances of He I are included. The shape of these resonances
  is given by a Fano profile - see Rumph, Bowyer, & Vennes 1994, AJ, 107, 2108.
  If these resonances are included then the input wavelength vector should have
  a fine (>~0.01 A) grid between 190 A and 210 A, since the resonances are very narrow.

### Output ###

* `tau`: Vector giving resulting optical depth, non-negative values.

### Example ###

One has a model EUV spectrum with wavelength, w (in Angstroms).
Find the EUV optical depth by 1e18 cm-2 of HI, with N(HeI)/N(HI) = N(HeII)/N(HI) = 0.05.

```jldoctest
julia> using AstroLib

julia> ismeuv.([670, 910], 1e19, 5e17, 5e17)
2-element Vector{Float64}:
 27.35393320556168
 62.683796028917286
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
ismeuv(wave::Real, hcol::Real, he1col::Real=hcol/10, he2col::Real=0,
       fano::Bool=false) =
           ismeuv(promote(float(wave), float(hcol), float(he1col),
                          float(he2col))..., fano)
