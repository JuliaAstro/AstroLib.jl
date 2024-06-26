# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function vactoair(wave_vac::AbstractFloat)
    if wave_vac >= 2000
        sigma2 = (1e4/wave_vac)^2  # Convert to wavenumber squared
        # Computer conversion factor.
        fact = 1 + 5.792105e-2/(238.0185 - sigma2) +
            1.67917e-3/(57.362 - sigma2)
        return wave_vac/fact # Convert wavelength
    else
        return wave_vac
    end
end

"""
    vactoair(wave_vacuum) -> wave_air

### Purpose ###

Converts vacuum wavelengths to air wavelengths.

### Explanation ###

Corrects for the index of refraction of air under standard conditions.
Wavelength values below ``2000 Å`` will not be altered.  Uses relation of
Ciddor (1996).

### Arguments ###

* `wave_vacuum`: vacuum wavelength in angstroms.  Wavelengths are corrected for the index of
  refraction of air under standard conditions.  Wavelength values below ``2000 Å`` will
  *not* be altered, take care within ``[1 Å, 2000 Å]``.

### Output ###

Air wavelength in angstroms.

### Method ###

Uses relation of Ciddor (1996), Applied Optics 35, 1566
(<http://adsabs.harvard.edu/abs/1996ApOpt..35.1566C>).

### Example ###

If the vacuum wavelength is `w = 2000`, then `vactoair(w)` yields an air
wavelength of `1999.353`.

```jldoctest
julia> using AstroLib

julia> vactoair(2000)
1999.3526230448367
```

### Notes ###

`airtovac` converts air wavelengths to vacuum wavelengths.

Code of this function is based on IDL Astronomy User's Library.
"""
vactoair(wave_vac::Real) = vactoair(float(wave_vac))
