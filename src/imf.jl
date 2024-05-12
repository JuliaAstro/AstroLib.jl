# This file is a part of AstroLib.jl. License is MIT "Expat".

function imf(mass::AbstractVector{T}, expon::AbstractVector{T},
             mass_range::AbstractVector{T}) where {T<:AbstractFloat}
    ne_comp = length(expon)
    if length(mass_range) != ne_comp + 1
        throw(ArgumentError(
            "Length of array mass_range is not one more than that of expon"))
    end
    integ = Vector{T}(undef, ne_comp)
    joint = ones(T, ne_comp)
    for i = 1:ne_comp
        if expon[i] != -1
            integ[i] = (mass_range[i+1]^(1 + expon[i]) - mass_range[i]^(1 + expon[i]))/
                       (1 + expon[i])
        else
            integ[i] = log(mass_range[i+1]/mass_range[i])
        end
        if i != 1
            joint[i] = joint[i-1]*(mass_range[i]^(expon[i-1] - expon[i]))
        end
    end
    norm = joint ./ (dot(integ, joint))
    psi = fill!(similar(mass), 0)
    for i = 1:ne_comp
        test = findall(mass_range[i] .<  mass .< mass_range[i+1])
        if length(test) !=0
            psi[test] = norm[i].*(mass[test].^expon[i])
        end
    end
    return psi
end

"""
    imf(mass, expon, mass_range) -> psi

### Purpose ###

Compute an N-component power-law logarithmic initial mass function (IMF).

### Explanation ###

The function is normalized so that the total mass distribution equals
one solar mass.

### Arguments ###

* `mass`: mass in units of solar mass, vector.
* `expon`: power law exponent, vector. The number of values in expon equals
  the number of different power-law components in the IMF.
* `mass_range`: vector containing the mass upper and lower limits of the
  IMF and masses where the IMF exponent changes. The number of values in
  mass_range should be one more than in expon. The values in mass_range
  should be monotonically increasing and positive.

### Output ###

* `psi`: mass function, number of stars per unit logarithmic mass interval
  evaluated for supplied masses.

### Example ###

Show the number of stars per unit mass interval at 3 Msun for a Salpeter
(expon = -1.35) IMF, with a mass range from 0.1 MSun to 110 Msun.

```jldoctest
julia> using AstroLib

julia> imf([3], [-1.35], [0.1, 110]) / 3
1-element Array{Float64,1}:
 0.01294143518151214
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
imf(mass::AbstractVector{<:Real}, expon::AbstractVector{<:Real}, mass_range::AbstractVector{<:Real}) =
    imf(float(mass), float(expon), float(mass_range))
