# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

# Helper function to format seconds part.
function formatsec(sec::Real, prec::Integer, truncate::Bool)
    sec = truncate ? trunc(sec, digits=prec) : round(sec, digits=prec)
    # Seconds of right ascension should be always positive (the hours part
    # holds the sign), so we don't need to take the absolute values.
    sec_frac, sec_int = modf(sec)
    # Format the integer part left padding with zeros.
    sec_int_str = lpad(string(trunc(Integer, sec_int)), 2, "0")
    # Unless precision is 0, format the fractional part with the decimal
    # separator "." followed by seconds rounded to the precision required
    # and right padded with zeros.
    sec_frac_str = prec == 0 ? "" : rpad(string(round(sec_frac, digits=prec+1)),
                                         prec+2, "0")[2:end]
    return sec_string = string(sec_int_str, sec_frac_str)
end

"""
    adstring(ra::Real, dec::Real[, precision::Int=2, truncate::Bool=true]) -> string
    adstring([ra, dec]) -> string
    adstring(dec) -> string
    adstring([ra], [dec]) -> ["string1", "string2", ...]

### Purpose ###

Returns right ascension and declination as string(s) in sexagesimal format.

### Explanation ###

Takes right ascension and declination expressed in decimal format, converts them
to sexagesimal and return a formatted string.  The precision of right ascension
and declination can be specified.

### Arguments ###

Arguments of this function are:

* `ra`: right ascension in decimal degrees.  It is converted to hours before
  printing.
* `dec`: declination in decimal degrees.

The function can be called in different ways:

* Two numeric arguments: first is `ra`, the second is `dec`.
* An iterable (array, tuple) of two elements: `(ra, dec)`.
* One numeric argument: it is assumed only `dec` is provided.

Optional keywords affecting the output format are always available:

* `precision` (optional integer keyword): specifies the number of digits of
  declination seconds.  The number of digits for right ascension seconds is
  always assumed to be one more `precision`.  If the function is called with
  only `dec` as input, `precision` default to 1, in any other case defaults to 0.
* `truncate` (optional boolean keyword): if true, then the last displayed digit
  in the output is truncated in precision rather than rounded.  This option is
  useful if `adstring` is used to form an official IAU name (see
  http://vizier.u-strasbg.fr/Dic/iau-spec.htx) with coordinate specification.

### Output ###

The function returns one string.  The format of strings can be specified with `precision`
and `truncate` keywords, see above.

### Example ###

```jldoctest
julia> using AstroLib

julia> adstring(30.4, -1.23, truncate=true)
" 02 01 35.9  -01 13 48"

julia> adstring.([30.4, -15.63], [-1.23, 48.41], precision=1)
2-element Vector{String}:
 " 02 01 36.00  -01 13 48.0"
 " 22 57 28.80  +48 24 36.0"
```
"""
function adstring(ra::T, dec::T; precision::Int=0,
                  truncate::Bool=false) where {T<:AbstractFloat}
    # XXX: IDL implementation takes also real values for "precision" and
    # truncates it.  I think it's better to enforce an integer type and cure
    # only negative values.
    precision = precision < 0 ? 0 : precision
    if isnan(ra)
        # If "ra" is NaN, print only declination.
        dec_deg, dec_min, dec_sec = sixty(dec)
        ra_string = ""
    else
        ra_hr, ra_min, ra_sec, dec_deg, dec_min, dec_sec = radec(ra, dec)
        ra_sec_string = formatsec(ra_sec, precision + 1, truncate)
        ra_string = @sprintf("%03.2d %02d %s  ", Int(ra_hr),
                             Int(ra_min), ra_sec_string)
    end
    dec_sec_string = formatsec(dec_sec, precision, truncate)
    dec_string = (dec >= 0 ? "+" : "-") * @sprintf("%02.2d %02d %s", Int(abs(dec_deg)), Int(dec_min), dec_sec_string)
    return string(ra_string, dec_string)
end

adstring(ra::Real, dec::Real;
         precision::Int=0, truncate::Bool=false) =
             adstring(promote(float(ra), float(dec))...,
                      precision=precision, truncate=truncate)

adstring(radec; precision::Int=0, truncate::Bool=false) =
    adstring(radec..., precision=precision, truncate=truncate)

# When printing only declination, IDL implementation defaults "precision" to 1
# instead of 0.
adstring(dec::Real; precision::Int=1, truncate::Bool=false) =
    adstring(NaN, dec, precision=precision, truncate=truncate)
