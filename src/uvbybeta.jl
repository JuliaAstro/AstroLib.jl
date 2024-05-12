# This file is a part of AstroLib.jl. License is MIT "Expat".

function _uvbybeta(by::T, m1::T, c1::T, hbeta::T, eby_in::T,
                   n::Integer) where {T<:AbstractFloat}
    # Rm1 = -0.33 & Rc1 = 0.19 & Rub = 1.53
    if n < 1 || n > 8
        throw(DomainError(
            "Input should be an integer in the range 1:8, giving approximate
            stellar classification"))
    end
    ub =  c1 + 2 * (m1 + by)
    # For group 1, beta is a luminosity indicator, c0 is a temperature indicator.
    # (u-b) is also a suitable temperature indicator.
    if n == 1
        # For dereddening, linear relation used between the intrinsic (b-y)
        # and (u-b) (Crawford 1978, AJ 83, 48)
        if isnan(eby_in)
           eby_in::T = (13.608 * by - ub + 1.467) / (12.078)
        end
        by0, m0, c0, ub0 = deredd(eby_in, by, m1, c1, ub)
        # When beta is not given, it is estimated using a cubic fit to the c0-beta
        # relation for luminosity class V given in Crawford (1978).
        if isnan(hbeta)
            hbeta::T = @evalpoly c0 2.61033 0.132557 0.161463 -0.027352
        end
        # Calculation of the absolute magnitude by applying the calibration of
        # Balona & Shobbrock (1974, MNRAS 211, 375)
        g = log10(hbeta - 2.515) - 1.6 * log10(c0 + 0.322)
        mv = 3.4994 + 7.2026 * log10(hbeta - 2.515) + @evalpoly(g, 0, -2.3192, 0, 2.9375)
        te = 5040 / (0.2917 * c0 + 0.2)
        # The ZAMS value of m0 is calculated from a fit to the data of Crawford (1978),
        # modified by Hilditch, Hill & Barnes (1983, MNRAS 204, 241)
        delm0 = @evalpoly(c0, 0.07473, 0.109804, -0.139003, 0.0957758) - m0
    elseif n == 2
        # For dereddening the linear relations between c0 and (u-b) determined from
        # Zhang (1983, AJ 88, 825) is used.
        if isnan(eby_in)
            eby_in = ((1.5 * c1 - ub + 0.035) / (1.5/(1.53/0.19) - 1)) / 1.53
        end
        by0, m0, c0, ub0 = deredd(eby_in, by, m1, c1, ub)

        if isnan(hbeta)
            hbeta = 2.542 + 0.037 * c0
        end
        delm0 = NaN
    elseif n == 3
        # For dereddening the linear relations between c0 and (u-b) determined from
        # Zhang (1983, AJ 88, 825) is used.
        if isnan(eby_in)
            eby_in = ((1.36 * c1 - ub + 0.004) / (1.36/(1.53/0.19) - 1)) / 1.53
        end
        by0, m0, c0, ub0 = deredd(eby_in, by, m1, c1, ub)
        # When beta is not given, it is derived from a fit of the c0-beta
        # relation of Zhang (1983).
        if isnan(hbeta)
            hbeta = 2.578 + 0.047 * c0
        end
        delm0 = NaN
    elseif n == 4
        # For dereddening the linear relations between c0 and (u-b) determined from
        # Zhang (1983, AJ 88, 825) is used.
        if isnan(eby_in)
            eby_in = ((1.32 * c1 - ub - 0.056) / (1.32/(1.53/0.19) - 1)) / 1.53
        end
        by0, m0, c0, ub0 = deredd(eby_in, by, m1, c1, ub)
        # When beta is not given, it is derived from a fit of the c0-beta
        # relation of Zhang (1983).
        if isnan(hbeta)
            hbeta = 2.59 + 0.066 * c0
        end
        delm0 = NaN
    elseif n == 5
        # For group 5, the hydrogen Balmer lines are at maximum; hence two new
        # parameters, a0 = f{(b-y),(u-b)} and r = f{beta,[c1]} are defined in
        # order to calculate absolute magnitude and metallicity.
        if isnan(eby_in)
            by0 = @evalpoly(m1 + 0.33 * by, -0.0235, -0.53921, 4.2608)
            while true
                bycorr = by0
                by0 = @evalpoly(m1 + 0.33 * (by - bycorr), 0.175709, -3.36225, 14.0881)

                if abs(bycorr - by0) < 0.001
                    break
                end
            end
            eby_in = by - by0
        end
        by0, m0, c0, ub0 = deredd(eby_in, by, m1, c1, ub)

        if isnan(hbeta)
            hbeta = 2.7905 - 0.6105 * by + 0.5 * m0 + 0.0355 * c0
        end
        a0 = by0 + 0.18 * (ub0 - 1.36)
        # MV is calculated according to Stroemgren (1966, ARA&A 4, 433)
        # with corrections by Moon & Dworetsky (1984, Observatory 104, 273)
        mv = 1.5 + 6 * a0 - 17 * (0.35 * (c1 - 0.19 * by) - (hbeta - 2.565))
        te = 5040 / (0.7536 * a0 + 0.5282)
        delm0 = @evalpoly(by0, 0.1598, 0.86888, -3.95105) - m0
    elseif n == 6

        if isnan(hbeta)
            hbeta = 3.06 - 1.221 * by - 0.104 * c1
        end
        m1zams = @evalpoly(hbeta, -17.209, 12.26, -2.158)

        if hbeta <= 2.74
            c1zams = 3 * hbeta - 7.56
            mvzams = 22.14 - 7 * hbeta
        elseif 2.74 < hbeta <= 2.82
            c1zams = 2 * hbeta - 4.82
            mvzams = 11.16 - 3 * hbeta
        else
            c1zams = 2 * hbeta - 4.83
            mvzams = @evalpoly hbeta -696.41 497.2 -88.4
        end

        if isnan(eby_in)
            delm1 = m1zams - m1
            delc1 = c1 - c1zams

            if delm1 < 0
                by0 = 2.946 - hbeta - 0.1 * delc1 - 0.25 * delm1
            else
                by0 = 2.946 - hbeta - 0.1 * delc1
            end
            eby_in = by - by0
        end
        by0, m0, c0, ub0 = deredd(eby_in, by, m1, c1, ub)
        delm0 =  m1zams - m0
        mv = mvzams - 9 * (c0 - c1zams)
        te = 5040 / (0.771453 * by0 + 0.546544)
    elseif n == 7
        # For group 7 c1 is the luminosity indicator for a particular beta, while
        # beta {or (b-y)0} indicates temperature.

        # Where beta is not available iteration is necessary to evaluate a corrected
        # (b-y) from which beta is then estimated.
        if isnan(hbeta)
            byinit = by
            m1init = m1
            for i in 1:1000
                m1by = @evalpoly byinit 0.345 -1.32 2.5
                bycorr = byinit + (m1by - m1init) / 2

                if abs(bycorr - byinit) <= 0.0001
                    break
                end
                byinit = bycorr
                m1init = m1by
            end
            hbeta = @evalpoly byinit 2.96618 -1.32861 1.01425
        end
        # m1(ZAMS) and mv(ZAMS) are calculated according to Crawford (1975) with
        # corrections suggested by Hilditch, Hill & Barnes (1983, MNRAS 204, 241)
        # and Olson (1984, A&AS 57, 443).
        m1zams = @evalpoly hbeta 46.4167 -34.4538 6.41701
        mvzams = @evalpoly hbeta 324.482 -188.748 11.0494 5.48012
        # c1(ZAMS) calculated according to Crawford (1975)
        if hbeta <= 2.65
            c1zams = 2 * hbeta - 4.91
        else
            c1zams = @evalpoly hbeta 72.879 -56.9164 11.1555
        end

        if isnan(eby_in)
            dbeta = 2.72 - hbeta
            eby_in = by - (0.222 - 0.05 * (c1 - c1zams) +
                    (1.11 - (0.1 + 3.6 * (m1zams - m1))) * dbeta + 2.7 * dbeta^2)
        end
        by0, m0, c0, ub0 = deredd(eby_in, by, m1, c1, ub)
        delm0 = m1zams - m0
        mv = mvzams - (9 + 20 * dbeta) * (c0 - c1zams)
        te = 5040/(0.771453 * by0 + 0.546544)
    elseif n == 8
        # Dereddening is done using color-color relations derived from
        # Olson (1984, A&AS 57, 443)
        if by <= 0.65
            eby_in = (5.8651 * by - ub - 0.8975) / 4.3351
        elseif 0.65 < by < 0.79
            eby_in = (-0.7875 * by - c1 + 0.6585) / -0.9775
            by0 = by - eby_in

            if by0 < 0.65
                eby_in = (5.8651 * by - ub - 0.8975) / 4.3351
            end
        else
            eby_in = (0.5126 * by - c1 - 0.3645) / 0.3226
            by0 = by - eby_in

            if by0 < 0.79
                eby_in = (-0.7875 * by - c1 + 0.6585) / -0.9775
            end
            by0 = by - eby_in
        end
        by0, m0, c0, ub0 = deredd(eby_in, by, m1, c1, ub)
        # m1(ZAMS), c1(ZAMS), and MV(ZAMS) are calculated according to Olson (1984)
        m1zams = @evalpoly by0 7.18436 -49.43695 122.1875 -122.466 42.93678

        if by0 < 0.65
            c1zams = @evalpoly by0 3.78514 -21.278 42.7486 -28.7056
            mvzams = @evalpoly by0 -59.2095 432.156 -1101.257 1272.503 -552.48
        elseif 0.65 <= by0 < 0.79
            c1zams = @evalpoly by0 0.33657 0.116031 -0.631821
            mvzams = @evalpoly by0 3.4305 4.97911 1.37632
        else
            c1zams = @evalpoly by0 -0.37237 0.530426 -0.010028
            mvzams = @evalpoly by0 4.37507 3.92776 1.18298
        end
        delm0 = m1zams - m0
        # Teff and MV calibration of Olson (1984)
        if by0 < 0.505
            f = 10 - 80 * (by0 - 0.38)
            te = exp10(3.924 - 0.416 * by0)
        else
            f = zero(T)
            te = exp10(3.869 -0.341 * by0)
        end
        mv = mvzams + f * (c1zams - c0) + 3.2 * delm0 - 0.07
    end

    if 1 < n < 5
        # c0-beta relation for ZAMS stars according to Crawford (1978,AJ 83, 48),
        # modified by Hilditch, Hill & Barnes (1983, MNRAS 204, 241).
        b_betaza = @evalpoly(c0, 2.62745, 0.228638, -0.099623, 0.277363, -0.160402) - 2.5
        # MV(ZAMS) calculated according to Balona & Shobbrock (1984, MNRAS 211, 375)
        mvzams = @evalpoly b_betaza -9.563 77.18 -206.98 203.704
        # MV is calculated from the d(beta)-d(MV) relation of Zhang (1983)
        dbeta = b_betaza - hbeta + 2.5
        mv = mvzams - @evalpoly(dbeta, 0.08, 61, -121.6)
        # Estimate of Teff by coupling the relations of Boehm-Vitense
        # (1981, ARA&A 19, 295) and Zhang (1983)
        te = 5040 / (0.27346 + 0.35866 * ub0)
    end
    # Transformation according to the FV-(b-y)0 relation of Moon (1984, MNRAS 211, 21P)
    if by0 <= 0.335
        fv = @evalpoly by0 3.981 -1.092 3.731 -6.759
    else
        fv = 3.959 - 0.534 * by0
    end
    radius = exp10(2 * (4.236 - 0.1 * mv - fv))
    eby = eby_in
    return te, mv, eby, delm0, radius
end

"""
    uvbybeta(by, m1, c1, n[, hbeta=NaN, eby_in=NaN]) -> te, mv, eby, delm0, radius

### Purpose ###

Derive dereddened colors, metallicity, and Teff from Stromgren colors.

### Arguments ###

* `by`: Stromgren b-y color
* `m1`: Stromgren line-blanketing parameter
* `c1`: Stromgren Balmer discontinuity parameter
* `n`: Integer which can be any value between 1 to 8, giving approximate stellar
  classification.
  (1) B0 - A0, classes III - V, 2.59 < Hbeta < 2.88,-0.20 <   c0   < 1.00
  (2) B0 - A0, class   Ia     , 2.52 < Hbeta < 2.59,-0.15 <   c0   < 0.40
  (3) B0 - A0, class   Ib     , 2.56 < Hbeta < 2.61,-0.10 <   c0   < 0.50
  (4) B0 - A0, class   II     , 2.58 < Hbeta < 2.63,-0.10 <   c0   < 0.10
  (5) A0 - A3, classes III - V, 2.87 < Hbeta < 2.93,-0.01 < (b-y)o < 0.06
  (6) A3 - F0, classes III - V, 2.72 < Hbeta < 2.88, 0.05 < (b-y)o < 0.22
  (7) F1 - G2, classes III - V, 2.60 < Hbeta < 2.72, 0.22 < (b-y)o < 0.39
  (8) G2 - M2, classes  IV - V, 0.20 < m0    < 0.76, 0.39 < (b-y)o < 1.00
* `hbeta` (optional): H-beta line strength index. If it is not supplied, then by
  default its value will be `NaN` and the code will estimate a value based on by,
  m1,and c1. It is not used for stars in group 8.
* `eby_in` (optional): specifies the E(b-y) color to use. If not supplied, then by
  default its value will be `NaN` and E(b-y) will be estimated from the Stromgren colors.

### Output ###

* `te`: approximate effective temperature
* `mv`: absolute visible magnitude
* `eby`: color excess E(b-y)
* `delm0`: metallicity index, delta m0, may not be calculable for early B stars
  and so returns `NaN`.
* `radius`: stellar radius (R/R(solar))

### Example ###

Suppose 5 stars have the following Stromgren parameters

by = [-0.001 ,0.403, 0.244, 0.216, 0.394]
m1 = [0.105, -0.074, -0.053, 0.167, 0.186]
c1 = [0.647, 0.215, 0.051, 0.785, 0.362]
hbeta = [2.75, 2.552, 2.568, 2.743, 0]
nn = [1,2,3,7,8]

Determine the stellar parameters

```jldoctest
julia> using AstroLib

julia> by = [-0.001 ,0.403, 0.244, 0.216, 0.394];

julia> m1 = [0.105, -0.074, -0.053, 0.167, 0.186];

julia> c1 = [0.647, 0.215, 0.051, 0.785, 0.362];

julia> hbeta = [2.75, 2.552, 2.568, 2.743, 0];

julia> nn = [1,2,3,7,8];

julia> uvbybeta.(by, m1, c1, nn, hbeta)
5-element Array{NTuple{5,Float64},1}:
 (13057.535222326893, -0.27375469585031265, 0.04954396423248884, -0.008292894218734928, 2.7136529525371897)
 (14025.053834219656, -6.907050783073221, 0.4140562248995983, NaN, 73.50771722263974)
 (18423.76405400214, -5.935816553877892, 0.2828247876690783, NaN, 39.84106215808709)
 (7210.507090112837, 2.2180408083364167, 0.018404079180028038, 0.018750927360588615, 2.0459018065648165)
 (5755.671513413262, 3.9449408311022, -0.025062997393370458, 0.03241423718769865, 1.5339239690774464)
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
uvbybeta(by::Real, m1::Real, c1::Real, n::Integer, hbeta::Real=NaN, eby_in::Real=NaN) =
    _uvbybeta(promote(float(by), float(m1), float(c1), float(hbeta), float(eby_in))..., n)
