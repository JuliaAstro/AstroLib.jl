# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function airtovac(wave_air::AbstractFloat)
    if wave_air < 2000
        return wave_air
    end
    wave_vac = wave_air
    for iter in 1:2
        sigma2 = (1.0e4 / wave_vac)^2 # Convert to wavenumber squared
        # Computer conversion factor.
        fact = 1 + 5.792105e-2 / (238.0185 - sigma2) +
            1.67917e-3 / (57.362 - sigma2)
        wave_vac = wave_air * fact # Convert Wavelength
    end
    return wave_vac
end

"""
    airtovac(wave_air) -> wave_vacuum

### Purpose ###

Converts air wavelengths to vacuum wavelengths.

### Explanation ###

Wavelengths are corrected for the index of refraction of air under standard conditions.
Wavelength values below ``2000 Å`` will *not* be altered, take care within
``[1 Å, 2000 Å]``.  Uses relation of Ciddor (1996).

### Arguments ###

* `wave_air`: the wavelength in air.

### Output ###

Vacuum wavelength in angstroms.

### Method ###

Uses relation of Ciddor (1996), Applied Optics 62, 958.

### Example ###

If the air wavelength is `w = 6056.125` (a Krypton line), then `airtovac(w)`
yields a vacuum wavelength of `6057.8019`.

```jldoctest
julia> using AstroLib

julia> airtovac(6056.125)
6057.801930991426
```

### Notes ###

`vactoair` converts vacuum wavelengths to air wavelengths.

Code of this function is based on IDL Astronomy User's Library.
"""
airtovac(wave_air::Real) = airtovac(float(wave_air))
