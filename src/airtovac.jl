# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _airtovac{T<:AbstractFloat}(wave_air::T)
    if wave_air >= 2000.0
        wave_vac = wave_air
        for iter= 1:2
            sigma2 = (1e4/wave_vac)^2.0 # Convert to wavenumber squared
            # Computer conversion factor.
            fact = 1.0 + 5.792105e-2/(238.0185e0 - sigma2) +
                1.67917e-3/(57.362e0 - sigma2)
            wave_vac = wave_air*fact # Convert Wavelength
        end
        return wave_vac
    else
        return wave_air
    end
end

"""
    airtovac(wave_air) -> wave_vacuum

### Purpose ###

Converts air wavelengths to vacuum wavelengths.

### Explanation ###

Wavelengths are corrected for the index of refraction of air under standard
conditions.  Wavelength values below \$2000 Å\$ will not be altered.  Uses
relation of Ciddor (1996).

### Arguments ###

* `wave_air`: can be either a scalar or an array of numbers.  Wavelengths are
corrected for the index of refraction of air under standard conditions.
Wavelength values below \$2000 Å\$ will *not* be altered, take care within \$[1
Å, 2000 Å]\$.

### Output ###

Vacuum wavelength in angstroms, same number of elements as `wave_air`.

### Method ###

Uses relation of Ciddor (1996), Applied Optics 62, 958.

### Example ###

If the air wavelength is `w = 6056.125` (a Krypton line), then `airtovac(w)`
yields a vacuum wavelength of `6057.8019`.

### Notes ###

`vactoair` converts vacuum wavelengths to air wavelengths.

Code of this function is based on IDL Astronomy User's Library.
"""
airtovac(wave_air::Real) = _airtovac(float(wave_air))
