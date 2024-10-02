# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _frebin(image::AbstractArray{T}, nsout::S, nlout::S, total::Bool) where {T<:AbstractFloat,S<:Integer}

    # Determine the size of the input array
    ns = size(image, 1)
    nl = length(image)/ns

    # Determine if the new sizes are integral factors of the original sizes
    sbox = ns/nsout
    lbox = nl/nlout

    # Contraction by an integral amount
    if (sbox == round(Int, sbox)) && (lbox == round(Int, lbox)) && (ns % nsout == 0) && (nl % nlout == 0)
        array_shaped = reshape(image, (Int(sbox), nsout, Int(lbox), nlout))
        return dropdims(sum(array_shaped, dims=(1,3)), dims=(1,3)) .* (total ? 1.0 : 1 ./ (sbox.*lbox))
    end

    # Expansion by an integral amount
    if (nsout % ns == 0) && (nlout % nl == 0)
        xindex = (1:nsout) / (nsout/ns)
        if isone(nl)  # 1D case, linear interpolation
            interpfunc = extrapolate(interpolate(image, BSpline(Constant())), Flat())
            return interpfunc(xindex) * (total ? sbox : 1.)
        end
        yindex = (1:nlout) / (nlout/nl)
        interpfunc = extrapolate(interpolate(image, BSpline(Constant())), Flat())
        return [interpfunc(x, y) for x in xindex, y in yindex] .* (total ? sbox.*lbox : 1.)
    end

    ns1 = ns-1
    nl1 = nl-1

    # Do 1D case separately
    if isone(nl)
        result = zeros(eltype(image), nsout)
        for i ∈ 0:nsout-1
            rstart = i*sbox                # starting position for each box
            istart = floor(Int, rstart)
            rstop = rstart + sbox          # ending position for each box
            istop = Int(clamp(floor(rstop), 0, ns1))
            frac1 = rstart-istart
            frac2 = 1.0 - (rstop-istop)

            # add pixel values from istart to istop and subtract fractional pixel from istart to start and
            # fractional pixel from rstop to istop

            result[i+1] = sum(image[istart+1:istop+1]) - frac1*image[istart+1] - frac2*image[istop+1]
        end
        return result .* (total ? 1.0 : 1 ./ (sbox.*lbox))
    end

    # Now, do 2D case
    # First, bin second dimension
    temp = zeros(eltype(image), ns, nlout)
    # Loop on output image lines
    for i ∈ 0:nlout-1
        rstart = i*lbox                # starting position for each box 
        istart = floor(Int, rstart)
        rstop = rstart + lbox
        istop = Int(clamp(floor(rstop), 0, nl1))
        frac1 = rstart-istart
        frac2 = 1.0 - (rstop-istop)

        # add pixel values from istart to istop and subtract fractional pixel from istart to start and
        # fractional pixel from rstop to istop

        if istart == istop
            temp[:,i+1] .= (1 .- frac1 .- frac2).*image[:,istart+1]
        else
            temp[:,i+1] .= dropdims(sum(image[:,istart+1:istop+1], dims=2), dims=2) .- frac1.*image[:,istart+1] .- frac2.*image[:,istop+1]
        end
    end
    temp = temp'
    # Bin in first dimension
    result = zeros(eltype(image), nlout, nsout)
    # Loop on output image samples
    for i ∈ 0:nsout-1
        rstart = i*sbox                # starting position for each box
        istart = floor(Int, rstart)
        rstop = rstart + sbox          # ending position for each box
        istop = Int(clamp(floor(rstop), 0, ns1))
        frac1 = rstart-istart
        frac2 = 1.0 - (rstop-istop)

        # add pixel values from istart to istop and subtract fractional pixel from istart to start and
        # fractional pixel from rstop to istop

        if istart == istop
            result[:,i+1] .= (1 .- frac1 .- frac2).*temp[:,istart+1]
        else
            result[:,i+1] .= dropdims(sum(temp[:,istart+1:istop+1], dims=2), dims=2) .- frac1.*temp[:,istart+1] .- frac2.*temp[:,istop+1]
        end
    end
    return transpose(result) .* (total ? 1.0 : 1 ./ (sbox.*lbox))

end

"""
    frebin(image, nsout, nlout=1, total=false) -> rebinned_image

### Purpose ###

Shrink or expand the size of an array an arbitrary amount using interpolation

### Arguments ###

* `image`: the array representing the image to be rebinned.
* `nsout`: number of samples in the output image, numeric scalar.
* `nlout` (optional): number of lines in the output image, numeric scalar (default = 1).
* `total` (optional boolean keyword): if true, the output pixels will be the 
           sum of pixels within the appropriate box of the input image.  
           Otherwise they will be the average. Use of the `total` keyword 
           conserves total counts.

### Output ###

The resized image is returned as the function result.

### Example ###

Suppose one has an 800 x 800 image array, im, that must be expanded to
a size 850 x 900 while conserving the total counts. The pixel values are
the sum of the x and y coordinates in the original image:

```jldoctest
julia> using AstroLib

julia> image = [x+y for x in 1:800, y in 1:800];

julia> image1 = frebin(image, 850, 900, total=true);

julia> sum(image)
512640000

julia> sum(image1)
5.126400000000241e8
```

### Notes ###

If the input image sizes are a multiple of the output image sizes
then `frebin` is equivalent to summing over the grid multiples for 
compression, and simple pixel duplication on expansion.
 
If the number of output pixels are not integers, the output image
size will be truncated to an integer.  The platescale, however, will
reflect the non-integer number of pixels.  For example, if you want to
bin a 100 x 100 integer image such that each output pixel is 3.1
input pixels in each direction use:
    n = 100/3.1   # 32.2581
    image_out = frebin(image,n,n)
 
The output image will be 32 x 32 and a small portion at the trailing
edges of the input image will be ignored.

### History ###

Adapted from May 1998 STIS  version, written D. Lindler, ACC
Added /NOZERO, use INTERPOLATE instead of CONGRID, June 98 W. Landsman  
Fixed for nsout non-integral but a multiple of image size  Aug 98 D.Lindler
DJL, Oct 20, 1998, Modified to work for floating point image sizes when
    expanding the image. 
Improve speed by addressing arrays in memory order W.Landsman Dec/Jan 2001

Code of this function is based on IDL Astronomy User's Library.
"""
function frebin(image::AbstractArray{R}, nsout::Real, nlout::Real=1; total::Bool=false) where {R<:Real}
    _frebin(float(image), promote(floor(Int, nsout), floor(Int, nlout))..., total)
end
