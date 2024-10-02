# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _fshift(image::AbstractArray{T}, Δx::T, Δy::T) where {T<:AbstractFloat}

    # Separate shift into an integer and fractional shift
    intx = floor(Int, Δx)
    inty = floor(Int, Δy)
    fracx = Δx - intx
    fracy = Δy - inty
    if fracx < 0
        fracx += 1
        intx -= 1
    end
    if fracy < 0
        fracy += 1
        inty -= 1
    end

    # Shift by the integer portion
    s = circshift(image, (intx, inty))
    if iszero(fracx) && iszero(fracy)
        return s
    end

    # Use bilinear interpolation between four pixels
    return s .* ((1 .- fracx) .* (1 .- fracy)) .+ 
           circshift(s, (0,1)) .* ((1 .- fracx) .* fracy) .+
           circshift(s, (1,0)) .* (fracx .* (1 .- fracy)) .+
           circshift(s, (1,1)) .* fracx .* fracy
end

"""
    fshift(image, Δx, Δy) -> shifted_image 

### Purpose ###

Routine to shift an image by non-integer values

### Arguments ###

* `image`: 2D image to be shifted.
* `Δx`: shift in x
* `Δy`: shift in y

### Output ###

Shifted image is returned as the function results

### Example ###

Suppose we want to shift a 10x10 image by 0.5 pixels in the
x and y directions. The pixel values are the sum of the x and
y coordinates:

```jldoctest
julia> using AstroLib

julia> image = [x+y for x in 1:10, y in 1:10];

julia> fshift(image, 0.5, 0.5)
10×10 Matrix{Float64}:
 11.0   7.0   8.0   9.0  10.0  11.0  12.0  13.0  14.0  15.0
  7.0   3.0   4.0   5.0   6.0   7.0   8.0   9.0  10.0  11.0
  8.0   4.0   5.0   6.0   7.0   8.0   9.0  10.0  11.0  12.0
  9.0   5.0   6.0   7.0   8.0   9.0  10.0  11.0  12.0  13.0
 10.0   6.0   7.0   8.0   9.0  10.0  11.0  12.0  13.0  14.0
 11.0   7.0   8.0   9.0  10.0  11.0  12.0  13.0  14.0  15.0
 12.0   8.0   9.0  10.0  11.0  12.0  13.0  14.0  15.0  16.0
 13.0   9.0  10.0  11.0  12.0  13.0  14.0  15.0  16.0  17.0
 14.0  10.0  11.0  12.0  13.0  14.0  15.0  16.0  17.0  18.0
 15.0  11.0  12.0  13.0  14.0  15.0  16.0  17.0  18.0  19.0
```

### History ###

version 2  D. Lindler  May, 1992 - rewritten for IDL version 2
19-may-1992	JKF/ACC		- move to GHRS DAF.

Code of this function is based on IDL Astronomy User's Library.

"""
function fshift(image::AbstractArray{R}, Δx::Real, Δy::Real) where {R<:Real}
    _fshift(float(image), promote(float(Δx), float(Δy))...)
end