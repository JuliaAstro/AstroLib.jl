# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

@testset "adstring" begin
    @test @inferred(adstring((30.4, -1.23), truncate=true)) ==
        @inferred(adstring([30.4, -1.23], truncate=true)) ==
        " 02 01 35.9  -01 13 48"
    @test @inferred(adstring(15, -.5)) == " 01 00 00.0  -00 30 00"
    @test @inferred(adstring(19.19321, truncate=true)) == "+19 11 35.5"
    @test @inferred(adstring(ten(36,24,15.015), -ten(8,24,36.0428), precision=3)) ==
        " 02 25 37.0010  -08 24 36.043"
    @test adstring.([30.4, -15.63], [-1.23, 48.41], precision=2) ==
        [" 02 01 36.000  -01 13 48.00", " 22 57 28.800  +48 24 36.00"]
    @test adstring.([(-58, 724)]) == [" 20 08 00.0  +724 00 00"]
end

@testset "airtovac" begin
    @test airtovac.([1234 6056.125]) ≈ [1234.0 6057.801930991426]
    @test @inferred(airtovac(2100)) ≈ 2100.666421596007
end

@testset "aitoff" begin
    @test @inferred(aitoff([227.23, 130], [-8.890, -35])) ==
        ([-137.92196683723276, 115.17541338020645],
         [-11.772527357473054, -44.491889962090085])
    @test @inferred(aitoff([375], [2.437])) ==
        ([16.63760711611838],[2.712427279646118])
    @test @inferred(aitoff((227.23, -8.890))) ==
        (-137.92196683723276,-11.772527357473054)
end

@testset "altaz2hadec" begin
    @test @inferred(altaz2hadec(59.086111, 133.30806, 43.07833)) ==
        @inferred(altaz2hadec((59.086111, 133.30806), 43.07833)) ==
        (336.68286017949157, 19.182449588316555)
    @test @inferred(altaz2hadec([15, 25, 35], [25.12, 45.32, -20.3],
                                [-23.44, 45.0, 52.5])) ==
        ([324.9881067314537, 256.7468302330436, 132.4919217875949] ,
         [44.38225395397647, 48.542947077386664, 67.33061196497327])
end

@testset "baryvel" begin
    dvelh_o, dvelb_o = @inferred(baryvel(2000))
    @test dvelh_o ≈ [1.582939967296732e-7, -1.0743272343303577e-7, -4.852410351888098e-8]
    @test dvelb_o ≈ [1.583299045307303e-7, -1.0736465601278539e-7, -4.849238001189245e-8]
    dvelh_o, dvelb_o = @inferred(baryvel(AstroLib.J2000, 1950))
    @test dvelh_o ≈ [-29.85888093436655, -4.684571288755146, -2.0305672776315777]
    @test dvelb_o ≈ [-29.849737241231153, -4.696440255370953, -2.035884790519881]
    dvelh_o, dvelb_o = @inferred(baryvel(jdcnv(1987, 4, 10, 0)))
    @test dvelh_o ≈ [9.514730366039178, -25.83862824722249, -11.202290688985778]
    @test dvelb_o ≈ [9.51331763088569,  -25.84985540772467, -11.207031996297271]
end

@testset "bprecess" begin
    ra, dec = @inferred(bprecess([ten(13, 42, 12.74)*15], [ten(8, 23, 17.69)],
                                 reshape(100*[-15*0.0257, -0.090], 2, 1)))
    @test ra  ≈ [204.93552515632123]
    @test dec ≈ [8.641287183886163]
    ra, dec = @inferred(bprecess(82, 19))
    @test ra  ≈ 81.26467916346334
    @test dec ≈ 18.959495700195394
    ra, dec = @inferred(bprecess([57], [23], 2024))
    @test ra  ≈ [56.26105898810067]
    @test dec ≈ [22.84693298145991]
    ra, dec = @inferred(bprecess([57], [23], reshape([9, 86], 2, 1), parallax=[1],
                                 radvel=[4]))
    @test ra  ≈ [56.25988479854577]
    @test dec ≈ [22.83493370392355]
    ra, dec = @inferred(bprecess((-57, -23), 2024))
    @test ra  ≈  302.2593299643789
    @test dec ≈ -23.150089972802036
    ra, dec = bprecess((-57, -23), [9, 86], parallax=1, radvel=4) # Inferred Type Error
    @test ra  ≈  302.2580376402947
    @test dec ≈ -23.16208183899836
end

@testset "calz_unred" begin
    @test calz_unred.(collect(900:1000:9900), ones(Float64, 10), -0.1) ≈
        [1.0, 0.43189326452379095, 0.5203675483533704, 0.594996469192435,
         0.6569506252451913, 0.7080829505773865, 0.7502392743978797, 0.7861262388745882,
         0.8151258710444882,0.8390325371659836]
end

# The values used for the testset are from running the code. They are slightly
# different from the output of the co_aberration routine of IDL AstroLib, as
# the function here uses an updated method to find mean obliquity
@testset "co_aberration" begin
    d_ra, d_dec = @inferred(co_aberration(jdcnv(1987, 4, 10, 0), ten(2,46,11.331)*15,
                                          ten(49,20,54.54), 1))
    @test d_ra ≈ -18.692441865574867
    @test d_dec ≈ -9.070782150537646
    ao, bo =  @inferred(co_aberration([57555.0, -6.44311e5], [302.282, 69.5667],
                                      [37.1519, 20.6847]))
    @test ao[1] ≈ 21.673056337579048
    @test ao[2] ≈ 18.496516329468466
    @test bo[1] ≈ -6.773070772568567
    @test bo[2] ≈ 2.9205843718089457
end

# The values used for the testset are from running the code. They are slightly
# different from the output of the co_aberration routine of IDL AstroLib, as
# the function here uses an updated method to find mean obliquity
@testset "co_nutate" begin
    an,bn,cn,dn,en = @inferred(co_nutate([jdcnv(2028,11,13,4,56), jdcnv(2013, 4, 16)],
                                         [10, 160],[80,30]))
    @test an ≈ [0.001209279097382776, 0.002465026741191423   ]
    @test bn ≈ [0.0017471459185713911, -0.0018134211486149354]
    @test cn ≈ [0.40904016038217567, 0.4090340058477726      ]
    @test dn ≈ [14.8593894278967, 12.102640377483143         ]
    @test en ≈ [2.7038090372351267, -5.86229256359996        ]
    ra_out, dec_out, eps_out, d_psi_out, d_eps_out = @inferred(co_nutate(2.451545e6,
                                                                         325, 0))
    @test ra_out ≈ -0.0035484441576727477
    @test dec_out ≈ -0.00034017946720967174
    @test eps_out ≈ 0.4090646078966446
    @test d_psi_out ≈ -13.923152677481191
    @test d_eps_out ≈ -5.773909654153591
end

@testset "co_refract" begin
    @test @inferred(co_refract(0.8)) ≈ 0.3714184384944585
    @test co_refract.([5.86,20], 50, 568.967, 273, 0.15, to_observe=true) ≈
        [5.94252628176525 , 20.02627520026167]
    @test @inferred(co_refract(14, 15000)) ≈ 13.990329255193124
end

@testset "ct2lst" begin
    @test ct2lst.(-76.72, -4, [DateTime(2008, 7, 30, 15, 53)]) ≈ [11.356505172312609]
    @test ct2lst.(9, [jdcnv(2015, 11, 24, 12, 21)]) ≈ [17.159574059885927]
end

# Test daycnv with Gregorian Calendar in force.
@testset "daycnv" begin
    @test @inferred(daycnv(2440000.0)) == DateTime(1968, 05, 23, 12)
    # Test daycnv with Julian Calendar in force (same result as IDL AstroLib's
    # daycnv).
    @test @inferred(daycnv(2000000.0)) == DateTime(763, 09, 18, 12)
    @test @inferred(daycnv(0.0)) == DateTime(-4713, 11, 24, 12)
end

@testset "deredd" begin
    by0, m0, c0, ub0 = @inferred(deredd([0.5, -0.5], [0.2, 0.5], [1, 1], [1, 1],
                                        [0.1, 0.3]))
    @test by0 ≈ [-0.3,0.5]
    @test m0  ≈ [1.165,1.0]
    @test c0  ≈ [0.905,1.0]
    @test ub0 ≈ [-0.665,0.3]
end

@testset "eci2geo" begin
    lat, long, alt = @inferred(eci2geo([0], [0], [0], [2452343]))
    @test lat  ≈ [0]
    @test long ≈ [12.992783145436988]
    @test alt  ≈ [-6378.137]
    lat, long, alt = @inferred(eci2geo((6978.137, 0, 0), jdcnv("2015-06-30T14:03:12.857")))
    @test lat  ≈ 0
    @test long ≈ 230.87301833205856
    @test alt  ≈ 600
    # Test `eci2geo' is the inverse of `geo2eci'
    jd = @inferred(get_juldate())
    lat, long, alt = @inferred(eci2geo(geo2eci(10, 10, 10, jd), jd))
    @test lat  ≈ 10
    @test long ≈ 10
    @test alt  ≈ 10
end

# The values used for the testset are from running the code. However they have been
# correlated with the output from eq2hor routine of IDL AstroLib, with
# differences only in the least significant digits.
@testset "eq2hor" begin
    alt_o, az_o, ha_o = eq2hor(259.20238631600944, 49.674907472176095,
                               AstroLib.J2000, "pbo", B1950=true)
    @test alt_o ≈ 43.804935603297004
    @test az_o ≈ 56.74141416977815
    @test ha_o ≈ 291.2750910119525
    alt_o, az_o, ha_o = eq2hor(142.2933457820434, -34.218006262991786, 2e6, 54.435, -34.78,
                               1000.34, ws=true, B1950=true, precession = false,
                               nutate=false, aberration=false,
                               refract=false, pressure = 500.345, temperature = 293.343)
    @test alt_o ≈ 1.3449999999999924
    @test az_o ≈ 359.43
    @test ha_o ≈ 359.3108663499664
    alt_o, az_o, ha_o = eq2hor(3.3222617779538037, 15.190516725395284, 2466879.7083333,
                               "kpno", pressure = 711, temperature = 273)
    @test alt_o ≈ 37.91138916818937
    @test az_o ≈ 264.918333213257
    @test ha_o ≈ 54.61193155973385
    alt_o, az_o, ha_o = @inferred(AstroLib._eq2hor(259.52076321839485, 49.62352289872951,
                                                   Float64(AstroLib.J2000), 43.0783,
                                                   -89.865, 0.0, NaN, NaN, false, false,
                                                   true, true, true, true))
    @test alt_o ≈ 43.687900264047116
    @test az_o ≈ 56.68399934960606
    @test ha_o ≈ 291.0817909922114
    alt_o, az_o = eq2hor(hor2eq(25, 55, 2.05e6, "pbo")[1:2]..., 2.05e6, "pbo")
    @test alt_o ≈ 24.99993224731665
    @test az_o ≈ 54.99993893556545
end

@testset "eqpole" begin
    x, y = @inferred(eqpole([100], [35], southpole=true))
    @test x ≈ [-111.18287262822456]
    @test y ≈ [ -19.604540237028665]
    x, y = @inferred(eqpole([80], [19]))
    @test x ≈ [72.78853915267848]
    @test y ≈ [12.83458333897169]
end

# The values used for the testset are from running the code. However they have been
# correlated with the output from euler routine of IDL AstroLib, with
# differences only in the least significant digits.
@testset "euler" begin
    glong, glat = @inferred(euler(299.590315, 35.201604, 1))
    @test glong ≈ 71.33498957116959
    @test glat ≈ 3.0668335310640984
    ra, dec = @inferred(euler((71.33498957116959, 3.0668335310640984), 2))
    @test ra ≈ 299.590315
    @test dec ≈ 35.201604
    elong, elat = @inferred(euler(3.141592653589793, 0.6143838917832061, 3,
                                  FK4 = true, radians=true))
    @test elong ≈ 2.8679433080257506
    @test elat ≈ 0.557258307291505
    ra, dec = @inferred(euler((2.8679433080257506, 0.557258307291505), 4,
                               FK4 = true, radians=true))
    @test ra ≈ 3.141592653589793
    @test dec ≈ 0.6143838917832061
    ecl, gal = euler.(30.45, 76.54, [5,6])
    @test ecl[1] ≈ 103.50477919192522
    @test ecl[2] ≈ 18.01965967759107
    @test gal[1] ≈ 194.96100731553986
    @test gal[2] ≈ 34.46136801388695
    @test @inferred(euler(183/pi, pi/180, 2, FK4=false, radians=true)) ==
        (5.682517110086799, 0.947078051715398)
    glong, glat = @inferred(euler([0.45, 130], [16.28, 53.65], 5))
    @test glong ≈ [96.9525940157568, 138.09922696730337]
    @test glat ≈ [-43.90672396295434, 46.95527026543361]
    @test_throws DomainError @inferred(euler((45,45), 7))
end

@testset "flux2mag" begin
    @test flux2mag.([1.5e-12, 8.7e-15, 4.4e-10]) ≈
        [8.459771852360795, 14.051201868453454, 2.291368308784527]
    @test @inferred(flux2mag(1)) ≈ -21.1
    @test @inferred(flux2mag(5.2e-15)) ≈ 14.609991640913002
    @test @inferred(flux2mag(5.2e-15, 15)) ≈ 20.709991640913003
    @test flux2mag(5.2e-15, ABwave=15) ≈ 27.423535345634598 # Inferred Type Error
end

@testset "gal_uvw" begin
    u, v, w = @inferred(gal_uvw([ten(1,9,42.3)*15], [ten(61,32,49.5)], [627.89], [77.84],
                                [-321.4], [1e3/129], lsr=true))
    @test u ≈ [118.2110474553902]
    @test v ≈ [-466.4828898385057]
    @test w ≈ [88.16573278565097]
    u, v, w = @inferred(gal_uvw(1, 2, 3, 4, 5, 6))
    @test u ≈  4.0228405867158745
    @test v ≈  3.7912174342038227
    @test w ≈ -3.1700191400725464
end

@testset "geo2eci" begin
    x, y, z = @inferred(geo2eci([0], [0], [0], [2452343]))
    @test x ≈ [6214.846433007192]
    @test y ≈ [-1433.9858454345972]
    @test z ≈ [0.0]
    x, y, z = @inferred(geo2eci((0,0,0), jdcnv("2015-06-30T14:03:12.857")))
    @test x ≈ -4024.8671780315185
    @test y ≈ 4947.835465127513
    @test z ≈ 0.0
end

@testset "geo2geodetic" begin
    lat, long, alt = @inferred(geo2geodetic([90], [0], [0], "Jupiter"))
    @test lat  ≈   [90]
    @test long ≈    [0]
    @test alt  ≈ [4638]
    lat, long, alt = @inferred(geo2geodetic((90, 0, 0)))
    @test lat  ≈ 90
    @test long ≈  0
    @test alt  ≈ 21.38499999999931
    lat, long, alt = @inferred(geo2geodetic((43.16, -24.32, 3.87), 8724.32, 8619.19))
    @test lat  ≈  43.849399515234516
    @test long ≈ -24.32
    @test alt  ≈  53.53354478670836
    lat, long, alt = @inferred(geo2geodetic([43.16], [-24.32], [3.87], 8724.32, 8619.19))
    @test lat  ≈ [ 43.849399515234516]
    @test long ≈ [-24.32]
    @test alt  ≈ [ 53.53354478670836]
end

@testset "geo2mag" begin
    lat, long = @inferred(geo2mag(ten(35,0,42), ten(135,46,6), 2016))
    @test lat  ≈  36.86579228937769
    @test long ≈ -60.184060536651614
    lat, long = @inferred(geo2mag([15], [24], 2016))
    @test lat  ≈ [ 11.452100529696096]
    @test long ≈ [-169.86030510727102]
end

@testset "geodetic2geo" begin
    lat, long, alt = @inferred(geodetic2geo([90], [0], [0], "Jupiter"))
    @test lat  ≈ [90]
    @test long == [0]
    @test alt  ≈ [-4638]
    lat, long, alt = @inferred(geodetic2geo((90, 0, 0)))
    @test lat ≈ 90
    @test long == 0
    @test alt ≈ -21.38499999999931
    lat, long, alt = @inferred(geodetic2geo((43.16, -24.32, 3.87), 8724.32, 8619.19))
    @test lat ≈ 42.46772711708433
    @test long == -24.32
    @test alt ≈ -44.52902080669082
    lat, long, alt = @inferred(geodetic2geo([43.16], [-24.32], [3.87], 8724.32, 8619.19))
    @test lat ≈ [42.46772711708433]
    @test long == [-24.32]
    @test alt ≈ [-44.52902080669082]
    # Test geodetic2geo is the inverse of geo2geodetic, within a certain
    # tolerance.
    lat, long, alt = @inferred(geodetic2geo(geo2geodetic(67.2,13.4,1.2)))
    @test lat ≈ 67.2 atol = 1e-8
    @test long == 13.4
    @test alt ≈ 1.2 atol = 1e-9
end

# Test get_date with mixed keywords.
@testset "get_date" begin
    @test @inferred(get_date(DateTime(2001,09,25,14,56,14), old=true,timetag=true)) ==
        @inferred(get_date(2001,09,25,14,56,14, old=true,timetag=true)) ==
        @inferred(get_date("2001-09-25T14:56:14", old=true,timetag=true)) ==
        "25/09/2001T14:56:14"
    @test @inferred(get_date(DateTime(2001,09,25,14,56,14))) ==
        @inferred(get_date(2001,09,25,14,56,14)) ==
        @inferred(get_date("2001-09-25T14:56:14")) ==
        "2001-09-25"
    @test get_date.([DateTime(2024), Date(2016, 3, 14)]) ==
        get_date.([Date(2024), "2016-03-14"]) ==
        get_date.(["2024-01", DateTime(2016, 3, 14)]) == ["2024-01-01", "2016-03-14"]
end

@testset "gcirc" begin
    @test gcirc.(0, [0,1,2], [1,2,3], [2,3,4], [3,4,5]) ≈
        [1.222450611061632, 2.500353926443337, 1.5892569925227757]
    @test @inferred(gcirc(0,  120, -43,   175, +22))     ≈  1.590442261600714
    @test @inferred(gcirc(1, (120, -43),  175, +22))     ≈  415908.56615322345
    @test @inferred(gcirc(2,  120, -43,  (175, +22)))    ≈  296389.3666794745
    @test @inferred(gcirc(0, (120, -43), (175, +22)) )   ≈  1.590442261600714
    @test gcirc.(1, [120], [-43],  175, +22)             ≈ [415908.56615322345]
    @test gcirc.(2,  120, -43,  [175], [+22])            ≈ [296389.3666794745]
    @test_throws DomainError @inferred(gcirc(3, 0, 0, 0, 0))
end

@testset "hadec2altaz" begin
    alt1, az1 = @inferred(hadec2altaz([0], [11.978165], [ten(43,4,42)]))
    @test alt1 ≈ [58.89983166666667]
    @test az1  ≈ [180.0]
    @test @inferred(hadec2altaz((0, 11.978165), ten(43,4,42), ws=true)[2]) ≈ 0.0
    alt1, az1 = 50, 20
    alt2, az2 = @inferred(hadec2altaz(altaz2hadec(alt1, az1, 40), 40))
    @test alt1 ≈ alt2
    @test az1  ≈ az2
end

@testset "helio_jd" begin
    @test helio_jd(juldate(2016, 6, 15, 11, 40), ten(20, 9, 7.8)*15, ten(37, 9, 7)) ≈
        57554.98808289718
    @test @inferred(helio_jd(1000, 23, 67, B1950=true)) ≈ 999.9997659545342
    @test @inferred(helio_jd(2000, 12, 88, diff=true))  ≈ -167.24845957792076
end

@testset "helio_rv" begin
    @test @inferred(helio_rv(helio_jd(juldate(94, 10, 25, 17, 30), ten(04, 38, 16)*15,
                                      ten(20, 41, 05)), 46487.5303, 2.0563056, -6, 59.3)) ≈
        -62.965570109145034
    @test helio_rv.([0.1, 0.9], 0, 1, 0, 100, 0.6, 45) ≈
        [-45.64994926111004, 89.7820347358485]
end

# The values used for the testset are from running the code. However they have been
# correlated with the output from helio routine of IDL AstroLib, with
# differences only in the least significant digits (except for `hrad`` output of Mars)
@testset "helio" begin
    @test_throws DomainError @inferred(helio(jdcnv(2005,07,17,2,6,9), 10))
    hrad_out, hlong_out, hlat_out = @inferred(helio(jdcnv(2000,08,23,0), 2, true))
    @test hrad_out ≈ 0.7213758288364316
    @test hlong_out ≈ 3.462574978561256
    @test hlat_out ≈ 0.050393862449261535
    hrad_out, hlong_out, hlat_out = @inferred(helio([AstroLib.J2000], [7]))
    @test hrad_out[1] ≈ 19.921687573575788
    @test hlong_out[1] ≈ 316.4011812518626
    @test hlat_out[1] ≈ -0.6846115653974465
end

# The values used for the testset are from running the code. However they have been
# correlated with the output from hor2eq routine of IDL AstroLib, with
# differences only in the least significant digits.
@testset "hor2eq" begin
    ra_o, dec_o, ha_o = hor2eq(43.6879, 56.684, AstroLib.J2000, "pbo", B1950=true)
    @test ra_o ≈ 259.3943636338339
    @test dec_o ≈ 49.67396411401468
    @test ha_o ≈ 291.0833116221485
    ra_o, dec_o, ha_o = hor2eq(1.345, 359.43, 2e6, 54.435, -34.78, 1000.34, ws=true,
                               B1950=true, precession = false, nutate=false,
                               aberration=false, refract=false, pressure = 500.345,
                               temperature = 293.343)
    @test ra_o ≈ 142.2933457820434
    @test dec_o ≈ -34.218006262991786
    @test ha_o ≈ 359.3108663499664
    ra_o, dec_o, ha_o = hor2eq(ten(37,54,41), ten(264,55,06), 2466879.7083333,
                               "kpno", pressure = 711, temperature = 273)
    @test ra_o ≈ 3.3222617779538037
    @test dec_o ≈ 15.190516725395284
    @test ha_o ≈ 54.61193186104758
    ra_o, dec_o, ha_o = @inferred(AstroLib._hor2eq(43.6879, 56.684, Float64(AstroLib.J2000),
                                                   43.0783, -89.865, 0.0, NaN, NaN, false,
                                                   false, true, true, true, true))
    @test ra_o ≈ 259.52076321839485
    @test dec_o ≈ 49.62352289872951
    @test ha_o ≈ 291.0817908419628
    ra_o, dec_o = hor2eq(eq2hor(45, 60, 2e6, "pbo")[1:2]..., 2e6, "pbo")
    @test ra_o ≈ 45.0001735428487
    @test dec_o ≈ 59.999995143349956
end

# The values used for the testset are from running the code. However they have been
# correlated with the output from imf routine of IDL AstroLib, with
# differences only in the least significant digits.
@testset "imf" begin
    @test_throws ArgumentError @inferred(imf([5], [-6.75], [0.9]))
    @test @inferred(imf([0.1, 0.01], [-0.6, -1], [ 0.007, 1.8, 110])) ≈
        [0.49627714725007616, 1.9757149090208912]
    @test imf.([[3],[5]], [[-1.35], [-0.6, -1.7]], [[0.1, 100], [0.007, 1.8, 110]]) ≈
        [[0.038948937298846235], [0.027349915327755464]]
end

# The values used for the testset are from running the code. However they have been
# correlated with the output from ismeuv routine of IDL AstroLib, with
# differences only in the least significant digits.
@testset "ismeuv" begin
    @test @inferred(ismeuv(304, 1e20)) ≈ 58.30508020244554
    @test ismeuv.([50, 75, 343], 1e18) ≈
        [0.004486567212077276, 0.01601219540569213, 0.7811526665834057]
    @test ismeuv.([96, 41, 233], 1e18, 1e17, 1e17) ≈
        [0.04657033979688484, 0.0035293162863089807, 0.30662192095758833]
    @test ismeuv.([96.6, 41.056, 233.19], 1e18, 1e17, 1e17, true) ≈
        [0.04733922036264192, 0.003543553203857853, 0.3103333860793942]
    @test ismeuv.([480, 910], 1e19, 5e17, 5e17) ≈
        [14.273937721143753, 62.68379602891728]
    @test @inferred(ismeuv(4500, 1e18)) ≈ 0
end

@testset "jdcnv" begin
    @test @inferred(jdcnv(-4713, 11, 24, 12)) ≈ 0.0
    @test @inferred(jdcnv(763, 09, 18, 12)) == @inferred(jdcnv("763-09-18T12")) == 2000000
    @test (jd=1234567.89; @inferred(jdcnv(daycnv(jd))) == jd)
    @test jdcnv.([DateTime(2016, 07, 31), "1969-07-20"]) ==
        jdcnv.([Date(2016, 07, 31), DateTime(1969, 07, 20)]) ==
        jdcnv.(["2016-07-31", Date(1969, 07, 20)])
end

@testset "jprecess" begin
    ra, dec = @inferred(jprecess([ten(13, 39, 44.526)*15], [ten(8, 38, 28.63)],
                                 reshape(100*[-15*0.0259, -0.093], 2, 1)))
    @test ra  ≈ [205.5530845731372]
    @test dec ≈ [8.388247441904628]
    ra, dec = @inferred(jprecess(82, 19))
    @test ra  ≈ 82.73568745151148
    @test dec ≈ 19.036972917272056
    ra, dec = @inferred(jprecess([57], [23], 2024))
    @test ra  ≈ [57.74049975335702]
    @test dec ≈ [23.150053754297726]
    ra, dec = @inferred(jprecess([57], [23], reshape([9, 86], 2, 1), parallax=[1],
                        radvel=[4]))
    @test ra  ≈ [57.74180294549785]
    @test dec ≈ [23.16200582079095]
    ra, dec = @inferred(jprecess((-57, -23), 2024))
    @test ra  ≈ 303.73910971499015
    @test dec ≈ -22.846895476784482
    ra, dec = jprecess((-57, -23), [9, 86], parallax=1, radvel=4) # Inferred Type Error
    @test ra  ≈ 303.7402950607101
    @test dec ≈ -22.834931625610313
end

# Test juldate with Gregorian Calendar in force.  This also makes sure precision
# of the result is high enough.  Note that "juldate(dt::DateTime) =
# Dates.datetime2julian(dt)-2.4e6" would not be precise.
@testset "juldate" begin
    @test @inferred(juldate(DateTime(2016, 1, 1, 8))) ≈ (57388.5 + 1//3)
    # Test juldate with Julian Calendar in force, for different centuries.  This
    # also makes sure precision of the result is high enough.
    @test @inferred(juldate(1582, 10, 1, 20))    ≈ (-100843 + 1//3)
    @test @inferred(juldate("1000-01-01T20"))    ≈ (-313692 + 1//3)
    @test @inferred(juldate("100-10-25T20"))     ≈ (-642119 + 1//3)
    @test @inferred(juldate(-4713, 1, 1, 12))    ≈ -2.4e6
    @test @inferred(juldate(2016, 06, 30, 00, 05, 53, 120)) ≈
        jdcnv(2016, 06, 30, 00, 05, 53, 120) - 2.4e6
    # Test daycnv and juldate together, with Gregorian Calendar in force.  Note that
    # they are not expected to be one the inverse of the other during Julian
    # Calendar.
    @test (dt=DateTime(2016, 1, 1, 20, 45, 33, 457);
           @inferred(daycnv(juldate(dt) + 2.4e6)) == dt)
end

@testset "kepler_solver" begin
    for e in 0:0.1:1
        M = collect(0:1e-3:2pi)
        E = mod2pi.(kepler_solver.(M, e))
        # Make sure E is the solution of the Kepler's Equation with high precision.
        @test M ≈ E .- e .* sin.(E) rtol = 1e-15
    end
    @test kepler_solver.([pi/4, pi/6, 8pi/3], 0) ≈ [pi/4, pi/6, 2pi/3]
    @test @inferred(kepler_solver(0, 1)) == 0.0
    @test_throws DomainError @inferred(kepler_solver(pi, -0.5))
    @test_throws DomainError @inferred(kepler_solver(pi,  1.5))
end

@testset "lsf_rotate" begin
    vel, lsf = @inferred(lsf_rotate(3, 90))
    @test @inferred(length(vel)) == @inferred(length(lsf)) == 61
    vel, lsf = @inferred(lsf_rotate(5, 10))
    @test vel ≈ collect(-10.0:5.0:10.0)
    @test lsf ≈ [0.0, 0.556914447710896, 0.6933098861837907, 0.556914447710896, 0.0]
end

@testset "mag2flux" begin
    @test @inferred(mag2flux(4.83, 21.12))               ≈ 4.1686938347033296e-11
    @test flux2mag(mag2flux(15, ABwave=12.), ABwave=12)  ≈ 15.0 # Inferred Type Error
    @test @inferred(mag2flux(8.3))                       ≈ 1.7378008287493692e-12
    @test @inferred(mag2flux(8.3, 12))                   ≈ 7.58577575029182e-9
    @test mag2flux(8.3, ABwave=12)                       ≈ 3.6244115683017193e-7 # Inferred Type Error
end

@testset "mag2geo" begin
    lat, long = @inferred(mag2geo(90, 0, 2016))
    @test lat  ≈ 86.395
    @test long ≈ -166.29000000000002
    lat, long = @inferred(mag2geo([15], [24], 2016))
    @test lat  ≈ [11.702066965890157]
    @test long ≈ [-142.6357492442842]
    # Test geo2mag is approximately the inverse of mag2geo
    lat, long = @inferred(geo2mag(mag2geo(12.34, 56.78, 2016)..., 2016))
    @test lat  ≈ 12.34
    @test long ≈ 56.78
end

@testset "mean_obliquity" begin
    @test @inferred(mean_obliquity(AstroLib.J2000)) ≈ 0.4090926006005829
    @test mean_obliquity.(jdcnv.([DateTime(1916, 09, 22, 03, 39),
                             DateTime(2063, 10, 13, 09)])) ≈
        [0.4092816887615259, 0.40894777540460037]
end

@testset "month_cnv" begin
    @test month_cnv.([" januavv  ", "SEPPES ", " aUgUsT", "la"]) == [1, 9, 8, -1]
    @test month_cnv.([2, 12, 6], short=true, low=true) == ["feb", "dec", "jun"]
    @test @inferred(month_cnv(5, up=true)) == "MAY"
    @test (list=[1, 2, 3]; month_cnv.(month_cnv.(list)) == list)
    @test (list=["July", "March", "November"]; month_cnv.(month_cnv.(list)) == list)
end

@testset "moonpos" begin
    ra, dec, dis, lng, lat = @inferred(moonpos(jdcnv(1992, 4, 12)))
    @test ra  ≈ 134.68846854844108
    @test dec ≈ 13.768366630560255
    @test dis ≈ 368409.68481612665
    @test lng ≈ 133.16726428105378
    @test lat ≈ -3.2291264192144356
    ra, dec, dis, lng, lat = @inferred(moonpos([2457521], radians=true))
    @test ra  ≈ [2.2587950290926178]
    @test dec ≈ [0.26183388011392217]
    @test dis ≈ [385634.68772395694]
    @test lng ≈ [2.232459255739553]
    @test lat ≈ [-0.059294466326164315]
end

@testset "mphase" begin
    @test mphase.([2457520, 2457530, 2457650]) ≈
        [0.2781695910737857, 0.9969808583803166, 0.9580708477591693]
end

@testset "nutate" begin
    long, obl = @inferred(nutate(jdcnv(1987, 4, 10)))
    @test long ≈ -3.787931077110755
    @test obl  ≈  9.442520698644401
    long, obl = @inferred(nutate(2457521))
    @test long ≈ -4.401443629818089
    @test obl  ≈ -9.26823431959121
    long, obl = @inferred(nutate([2457000, 2458000]))
    @test long ≈ [ 4.327189321653877, -9.686089990639474]
    @test obl  ≈ [-9.507794266102866, -6.970768250588256]
end

@testset "paczynski" begin
    @test @inferred(paczynski(-1e-10)) ≈  -1e10
    @test @inferred(paczynski(1e-1))   ≈  10.037461005722337
    @test @inferred(paczynski(-1))     ≈  -1.3416407864998738
    @test @inferred(paczynski(10))     ≈   1.0001922892047386
    @test @inferred(paczynski(-1e10))  ≈  -1
end

@testset "planck_freq" begin
    @test @inferred(planck_freq(2000, 5000)) ≈ 6.1447146126144004e-30
end

@testset "planck_wave" begin
    @test @inferred(planck_wave(2000, 5000)) ≈ 8.127064833530511e-24
end

# The values used for the testset are from running the code. However they have been
# correlated with the output from planet_coord routine of IDL AstroLib, with
# differences only in the least significant digits
@testset "planet_coords" begin
    @test_throws DomainError @inferred(planet_coords(DateTime(2013, 07, 22,
                                                              03, 19, 06),0))
    ra_out, dec_out = @inferred(planet_coords([AstroLib.J2000, 2.45e6], [2,8]))
    @test ra_out[1] ≈ 239.8965221579066
    @test ra_out[2] ≈ 294.55483837772476
    @test dec_out[1] ≈ -18.450868549676304
    @test dec_out[2] ≈ -20.992319312642262
    ra_out, dec_out = @inferred(planet_coords(2.45e6, 9))
    @test ra_out ≈ 238.3131013864547
    @test dec_out ≈ -6.964788781133789
    @test @inferred(planet_coords(juldate(1), 3) == (0, 0))
end

@testset "polrec" begin
    x, y = @inferred(polrec([1, 2, 3], [pi, pi/2.0, pi/4.0]))
    @test x ≈ [-1.0, 0.0, 1.5*sqrt(2.0)]
    @test y ≈ [ 0.0, 2.0, 1.5*sqrt(2.0)]
    x, y = @inferred(polrec((2, 135), degrees=true))
    @test x ≈ -sqrt(2)
    @test y ≈  sqrt(2)
end

@testset "posang" begin
    @test @inferred(posang(1, ten(13, 25, 13.5), ten(54, 59, 17),
                           ten(13, 23, 55.5), ten(54, 55, 31))) ≈ -108.46011246802047
    @test posang.(0, [0,1,2], [1,2,3], [2,3,4], [3,4,5]) ≈
        [1.27896824717634, 1.6840484573313608, 0.2609280020139511]
    @test @inferred(posang(0,  120, -43,   175, +22))     ≈ -1.5842896165356724
    @test @inferred(posang(1, (120, -43),  175, +22))     ≈ 82.97831348792039
    @test @inferred(posang(2,  120, -43,  (175, +22)))    ≈ 50.02816530382374
    @test @inferred(posang(0, (120, -43), (175, +22)))    ≈ -1.5842896165356724
    @test posang.(1, [120], [-43],  175, +22)  ≈ [82.97831348792039]
    @test posang.(2,  120, -43,  [175], [+22]) ≈ [50.02816530382374]
    @test_throws DomainError @inferred(posang(3, 0, 0, 0, 0))
end

@testset "precess" begin
    ra1, dec1 = @inferred(precess((ten(2,31,46.3)*15, ten(89,15,50.6)), 2000, 1985))
    @test ra1  ≈ 34.09470328718033
    @test dec1 ≈ 89.19647174928589
    ra2, dec2 = @inferred(precess([ten(21, 59, 33.053)*15], [ten(-56, 59, 33.053)],
                                  1950, 1975, FK4=true))
    @test ra2  ≈ [330.3144305418865]
    @test dec2 ≈ [-56.87186126487889]
end

# The values used for the testset are from running the code. However they have been
# correlated with the output from precess_cd routine of IDL AstroLib, with
# differences only in the least significant digits.
@testset "precess_cd" begin
    @test @inferred(precess_cd([30 60; 60 90], 1950, 2000, [13, 8], [43, 23])) ≈
        [30.919029003435927 62.343060521017435;
         61.93905850970097 93.56509103294071  ]
    @test @inferred(precess_cd([30 60; 60 90], 2000, 1950, [13, 8], [43, 23])) ≈
        [30.919029003435927 62.343060521017435;
         61.93905850970097 93.56509103294071  ]
    @test @inferred(precess_cd([12.45 56.7; 66 89], 2000, 1985, [67.4589455, 0.345345],
                               [37.94291666666666, 89.26405555555556])) ≈
        [963.4252080520984 4387.890452343343  ;
         5107.504395433958 6887.55936949333   ]
    @test @inferred(precess_cd([30.0 28.967; 60.45 90.65], 2000, 1975, [13, 10.658],
                               [35.54, 67], true)) ≈
        [64.78429186351575 62.637156996728194 ;
         130.49379143419722 195.9699513801844 ]
end

@testset "precess_xyz" begin
    x1, y1, z1 = @inferred(precess_xyz((1.2, 2.3, 1.7), 2000, 2050))
    @test x1 ≈ 1.165933061423247
    @test y1 ≈ 2.313228746401996
    @test z1 ≈ 1.7057470102860104
    x2, y2, z2 = @inferred(precess_xyz([0.7, -2.4], [3.3, 6.6], [0, 4], 2000, 2016))
    @test x2 ≈ [0.688187142071843,   -2.429815562246262]
    @test y2 ≈ [3.3024835038223532,   6.591359330834213]
    @test z2 ≈ [0.001079105285993004, 3.9962455511755794]
end

@testset "premat" begin
    @test @inferred(premat(1967, 1982, FK4=true)) ≈
        [0.9999933170034135    -0.0033529069683496567 -0.0014573823699636742;
         0.00335290696825777    0.9999943789886484    -2.443304965138481e-6 ;
         0.0014573823701750721 -2.4431788671274868e-6  0.9999989380147651   ]
    @test @inferred(premat(1995, 2003)) ≈
        [ 0.9999980977132219    -0.0017889257711428855 -0.0007773766929507687;
          0.0017889257711354528  0.9999983998707707    -6.953448226403318e-7 ;
          0.0007773766929678732 -6.953257000046125e-7   0.9999996978424512   ]
end

@testset "radec" begin
    @test @inferred(radec(15.90, -0.85)) == (1.0, 3.0, 36.0, -0.0, 51.0, 0.0)
    @test @inferred(radec(-0.85,15.9)) == (23.0,56.0,36.0,15.0,54.0,0.0)
    @test @inferred(radec(-20,4,hours=true)) == (4.0,0.0,0.0,4.0,0.0,0.0)
    @test @inferred(radec([15.90, -0.85], [-0.85,15.9])) ==
        ([1.0, 23.0], [3.0, 56.0], [36.0, 36.0],
         [-0.0, 15.0], [51.0, 54.0], [0.0, 0.0])
end

@testset "recpol" begin
    r = a = zeros(Float64, 3)
    r, a = @inferred(recpol([0, sqrt(2.0), 2.0*sqrt(3.0)], [0, sqrt(2.0), 2.0]))
    @test r ≈ [0.0,  2.0,  4.0]
    @test a ≈ [0.0, pi/4.0, pi/6.0]
    r, a = @inferred(recpol(1, 1))
    @test r ≈ sqrt(2.0)
    @test a ≈ pi/4.0
    # Test polrec is the inverse of recpol
    xi, yi, x, y, = 6.3, -2.7, 0.0, 0.0
    x, y = @inferred(polrec(recpol((xi, yi), degrees=true), degrees=true))
    @test x ≈ xi
    @test y ≈ yi
end

@testset "rhotheta" begin
    ρ, θ = @inferred(rhotheta(41.623, 1934.008, 0.2763, 0.907, 59.025, 23.717, 219.907,
                              1980))
    @test ρ ≈ 0.41101776646245836
    @test θ ≈ 318.4242564860495
end

# Test also it's the inverse of ten.
@testset "sixty" begin
    @test @inferred(sixty(-51.36))                    ≈ [-51.0, 21.0, 36.0]
    @test @inferred(ten(sixty(-0.10934835545824395))) ≈ -0.10934835545824395
    @test @inferred(sixty(1))                         ≈ [1.0, 0.0, 0.0]
end

@testset "sphdist" begin
    @test sphdist.([0,1,2], [1,2,3], [2,3,4], [3,4,5]) ≈
        [1.222450611061632, 2.500353926443337, 1.5892569925227762]
    @test @inferred(sphdist(120, -43, 175, +22))      ≈  1.5904422616007134
    @test sphdist.([120], [-43], 175, +22) ≈ [1.5904422616007134]
    @test sphdist.(120, -43, [175], [+22]) ≈ [1.5904422616007134]
end

@testset "sunpos" begin
    ra, dec, lon, obl = @inferred(sunpos(jdcnv(1982, 5, 1)))
    @test ra  ≈ 37.88589057369026
    @test dec ≈ 14.909699471099517
    @test lon ≈ 40.31067053890748
    @test obl ≈ 23.440840980112657
    ra, dec, lon, obl = @inferred(sunpos(jdcnv.([DateTime(2016, 5, 10)]), radians=true))
    @test ra  ≈ [0.8259691339090751]
    @test dec ≈ [0.3085047454107549]
    @test lon ≈ [0.8687853454154388]
    @test obl ≈ [0.40901175207670365]
    ra, dec, lon, obl = @inferred(sunpos([2457531]))
    @test ra  ≈ [59.71655864208797]
    @test dec ≈ [20.52127006818727]
    @test lon ≈ [61.824436793991545]
    @test obl ≈ [23.434648653514724]
end

# Make sure string and numerical inputs are consistent (IDL implementation of "ten" is not).
@testset "ten" begin
    @test @inferred(ten(0, -23, 34)) == @inferred(ten((0, -23, 34))) ==
        @inferred(ten([0, -23, 34])) == ten(" : 0 -23 :: 34") == -0.37388888888888894
    @test @inferred(ten(-0.0, 60)) == @inferred(ten((-0.0, 60))) ==
        @inferred(ten([-0.0, 60])) == @inferred(ten("-0.0 60")) == -1.0
    @test @inferred(ten(-5, -60, -3600)) == @inferred(ten((-5, -60, -3600))) ==
        @inferred(ten([-5, -60, -3600])) == @inferred(ten("  -5: :-60: -3600")) == -3.0
    @test ten("") == 0.0
    @test ten.([0, -0.0, -5], [-23, 60, -60], [34, 0, -3600]) ==
        ten.([(0, -23,34), ":-0.0:60", (-5, -60, -3600)]) ==
        ten.(["0   -23 :: 34", (-0.0, 60), " -5:-60: -3600"]) ==
        [-0.37388888888888894, -1.0, -3.0]
    @test ten.([12.0, -0.0], [24, 30]) == ten.([" 12::24", " -0:30: "]) == [12.4, -0.5]
end

@testset "tic_one" begin
    min2, tic1 = @inferred(tic_one(30.2345, 12.74, 10))
    @test min2 ≈ 30.333333333333332
    @test tic1 ≈ 7.554820000000081
    min2, tic1 = @inferred(tic_one(45, 50, 4, true))
    @test min2 ≈ 46.0
    @test tic1 ≈ 50.0
    min2, tic1 = @inferred(tic_one(pi\8, tics(90, 45, 1000, 10)...))
    @test min2 ≈ 2.5
    @test tic1 ≈ 1.0318357862412286
end

@testset "ticpos" begin
    @test ticpos.([16,8,4],[1024,512,256], [150,75,37.5]) ==
                  [(256.0, 4, "Degrees"), (128.0, 2, "Degrees"), (64.0, 1, "Degrees")]
    @test @inferred(ticpos(2, 512, 75)) == (128.0, 30, "Arc Minutes")
    @test @inferred(ticpos(1.5, 512, 75)) == (85.33333333333333, 15, "Arc Minutes")
    @test @inferred(ticpos(1.5, 512, 50)) == (56.888888888888886, 10, "Arc Minutes")
    @test @inferred(ticpos(1.6, 1024, 50)) == (53.333333333333336, 5, "Arc Minutes")
    @test @inferred(ticpos(0.2, 512, 75)) == (85.33333333333333, 2, "Arc Minutes")
    @test @inferred(ticpos(0.5, 512, 10)) == (17.066666666666666, 1, "Arc Minutes")
    @test @inferred(ticpos(0.1, 1024, 50)) == (85.33333333333333, 30, "Arc Seconds")
    @test @inferred(ticpos(0.08, 1024, 40)) == (53.333333333333336, 15, "Arc Seconds")
    @test @inferred(ticpos(0.025, 512, 50)) == (56.888888888888886, 10, "Arc Seconds")
    @test @inferred(ticpos(pi/100, 1024, 40)) == (45.27073936836133, 5, "Arc Seconds")
    @test @inferred(ticpos(0.06, 2048, 20)) == (18.962962962962965, 2, "Arc Seconds")
    @test @inferred(ticpos(0.016, 1024, 20)) == (17.77777777777778, 1, "Arc Seconds")
end

@testset "tics" begin
    @test @inferred(tics(30, 90, 30, 1)) == (3.8666666666666667, 480)
    @test @inferred(tics(30, 90, 3, 3, true)) == (4.0, 240)
    @test @inferred(tics(30, 70, 3, 1, true)) == (0.75, 60)
    @test tics.([30,50],[70,60], [6,12], [3,0.5], [true, false]) ==
                [(3.75,120), (0.55,30)]
    @test @inferred(tics(45, 55, 30, 0.5)) == (0.725, 15)
    @test @inferred(tics(45, 60, 10, 0.1)) == (0.1, 10)
    @test @inferred(tics(55, 60, 100.0, 1/2)) == (0.66, 2)
    @test @inferred(tics(25, 30, 50, 2, true)) == (2.45, 1)
    @test @inferred(tics(20, 80, 600, 0.03)) == (0.04159722222222222, 0.25)
    @test @inferred(tics(25, 75, 500, 0.02)) == (0.02772222222222222, 0.16666666666666666)
    @test @inferred(tics(10, 12, 25, 0.01)) == (0.016666666666666666, 0.08333333333333333)
    @test @inferred(tics(20, 80, 6000, 0.03)) ==
        (0.055546296296296295, 0.03333333333333333)
    @test @inferred(tics(30, 60, 200, 0.02, true)) ==
        (0.02763888888888889, 0.016666666666666666)
    @test @inferred(tics(60, 70, 125, 0.001)) ==
        (0.0017222222222222222, 0.008333333333333333)
    @test @inferred(tics(10, 12, 25, 0.01, true)) == (0.01, 0.0033333333333333335)
    @test @inferred(tics(130, 180, 1000, 0.0004)) == (0.000555, 0.0016666666666666668)
    @test @inferred(tics(60, 70, 5500//2, 0.003)) ==
        (0.003818055555555556, 0.0008333333333333334)
    @test @inferred(tics(30, 150, 4000, 0.002, true)) ==
        (0.0027770833333333337, 0.0003333333333333333)
    @test @inferred(tics(9.5, 14.5, 5000, 0.002)) ==
        (0.002777222222222222, 0.00016666666666666666)
    @test @inferred(tics(90, 45, 1000, 10)) == (11.1, -30)
    ticsize, incr = @inferred(tics(30, 70, 50, 0.1))
    @test ticsize ≈ 0.10208333333333333
    @test incr ≈ 5
    ticsize, incr = @inferred(tics(pi/3, pi/2, 60.0, 7.5, true))
    @test ticsize ≈ 14.085212463632736
    @test incr ≈ 0.5
end

@testset "true_obliquity" begin
    @test @inferred(true_obliquity(AstroLib.J2000)) ≈ 0.4090646078966446
    @test true_obliquity.(jdcnv.([DateTime(2016, 08, 23, 03, 39, 06),
                             DateTime(763, 09, 18, 12)])) ≈
        [0.4090133706884892, 0.41188965892279295]
end

@testset "kepler_solver" begin
    @test @inferred(trueanom(8pi/3, 0.7))              ≈ 2.6657104039293764
    @test trueanom.([pi/4, pi/6, 8pi/3], 0)            ≈ [pi/4, pi/6, 2pi/3]
    @test @inferred(trueanom(3pi/2, 0.8))              ≈ -2.498091544796509
    @test @inferred(trueanom(0.1, 1))                  ≈ pi
    @test_throws DomainError @inferred(trueanom(pi, -0.5))
    @test_throws DomainError @inferred(trueanom(pi,  1.5))
end

# The values used for the testset are from running the code. However they have been
# correlated with the output from uvbybeta routine of IDL AstroLib, with
# differences only in the least significant digits.
@testset "uvbybeta" begin
    @test_throws DomainError @inferred uvbybeta(NaN, NaN, NaN, 9, NaN)
    te_o, mv_o, eby_o, delm_o, radius_o = @inferred uvbybeta(-0.009, 0.11, 0.68, 1)
    @test te_o ≈ 12719.770257555097
    @test mv_o ≈ -0.08913388245565224
    @test eby_o ≈ 0.038295081967213124
    @test delm_o ≈ -0.007787863429965028
    @test radius_o ≈ 2.549253655861009
    te_o, mv_o, eby_o, delm_o, radius_o = @inferred uvbybeta(0.415, -0.069, 0.308, 2)
    @test te_o ≈ 12280.215047866226
    @test mv_o ≈ -7.112880191374548
    @test eby_o ≈ 0.4040160642570281
    @test isnan(delm_o)
    @test radius_o ≈ 90.29892932880738
    te_o, mv_o, eby_o, delm_o, radius_o = @inferred uvbybeta(0.273, -0.051, 0.051, 3)
    @test te_o ≈ 18733.204098765906
    @test mv_o ≈ -5.453550292618578
    @test eby_o ≈ 0.33158225857187795
    @test isnan(delm_o)
    @test radius_o ≈ 27.825306168356892
    te_o, mv_o, eby_o, delm_o, radius_o = @inferred uvbybeta(0.29343, 0.09, 0.001, 4)
    @test te_o ≈ 25713.1878000067
    @test mv_o ≈ -4.816603735163776
    @test eby_o ≈ 0.6430112570356474
    @test isnan(delm_o)
    @test radius_o ≈ 0.16611917129253512
    te_o, mv_o, eby_o, delm_o, radius_o = @inferred uvbybeta(0.03, 0.21345, 0.3976, 5)
    @test te_o ≈ 10363.815662133818
    @test mv_o ≈ 4.411076600000009
    @test eby_o ≈ -0.038394585955939545
    @test delm_o ≈ -0.03113954500000002
    @test radius_o ≈ 0.4863555386334208
    te_o, mv_o, eby_o, delm_o, radius_o = @inferred uvbybeta(0.18547, -0.04356, 0.2576, 6)
    @test te_o ≈ 7308.315628104092
    @test mv_o ≈ 7.562860949999999
    @test eby_o ≈ -0.00736941600000007
    @test delm_o ≈ 0.2449243827527609
    @test radius_o ≈ 0.17062904371865295
    te_o, mv_o, eby_o, delm_o, radius_o = @inferred uvbybeta(0.18547, -0.04356, 0.2576, 6,
                                                             2.72001)
    @test te_o ≈ 7308.315628104092
    @test mv_o ≈ 6.1817999999999955
    @test eby_o ≈ -0.07476300000000027
    @test delm_o ≈ 0.21601800458420034
    @test radius_o ≈ 0.3223035928541916
    te_o, mv_o, eby_o, delm_o, radius_o = @inferred uvbybeta(0.18547, -0.04356, 0.2576, 6,
                                                             2.82)
    @test te_o ≈ 7334.8182052629045
    @test mv_o ≈ 7.767123299999994
    @test eby_o ≈ 0.003229999999999733
    @test delm_o ≈ 0.2454149000000012
    @test radius_o ≈ 0.154376398171302
    te_o, mv_o, eby_o, delm_o, radius_o = @inferred uvbybeta(0.18547, -0.04356, 0.2576, 6,
                                                             2.8201)
    @test te_o ≈ 7343.722680158638
    @test mv_o ≈ 7.682471615999958
    @test eby_o ≈ 0.0043099999999998695
    @test delm_o ≈ 0.24506736642000115
    @test radius_o ≈ 0.16018990575428574
    te_o, mv_o, eby_o, delm_o, radius_o = @inferred uvbybeta(0.1923, 0.2186, 0.2783, 6)
    @test te_o ≈ 7252.900460401587
    @test mv_o ≈ 7.219177499999998
    @test eby_o ≈ -0.011614630064603726
    @test delm_o ≈ -0.019005720258414266
    @test radius_o ≈ 0.2024586440816956
    te_o, mv_o, eby_o, delm_o, radius_o = @inferred uvbybeta(0.216, 0.167, 0.785, 7)
    @test te_o ≈ 7108.688707930596
    @test mv_o ≈ 1.538642957897997
    @test eby_o ≈ 0.005426551873147473
    @test delm_o ≈ 0.010495724945951107
    @test radius_o ≈ 2.8675336941552185
    te_o, mv_o, eby_o, delm_o, radius_o = @inferred uvbybeta(0.206, 0.162, 0.786, 7,
                                                             2.601)
    @test te_o ≈ 7144.241056060125
    @test mv_o ≈ -0.8920873665838025
    @test eby_o ≈ -0.12713573819611615
    @test delm_o ≈ 0.05272866901000031
    @test radius_o ≈ 8.706345749177931
    te_o, mv_o, eby_o, delm_o, radius_o = @inferred uvbybeta(0.394, 0.184, 0.382, 8)
    @test te_o ≈ 5755.671513413262
    @test mv_o ≈ 3.7737408311021916
    @test eby_o ≈ -0.028753800373693867
    @test delm_o ≈ 0.034414237187698926
    @test radius_o ≈ 1.6597544196694192
    te_o, mv_o, eby_o, delm_o, radius_o = @inferred uvbybeta(0.6501, 0.184, 0.382, 8)
    @test te_o ≈ 5452.352129194577
    @test mv_o ≈ 5.192225306337874
    @test eby_o ≈ 0.19958051948051947
    @test delm_o ≈ 0.032731030075137835
    @test radius_o ≈ 0.9924525316819005
    te_o, mv_o, eby_o, delm_o, radius_o = @inferred uvbybeta(0.7901, 0.184, 0.0145, 8)
    @test te_o ≈ 3977.185981725903
    @test mv_o ≈ 9.893069966897801
    @test eby_o ≈ -0.022297953964194334
    @test delm_o ≈ 0.5456850156962505
    @test radius_o ≈ 0.26255047396110764
    te_o, mv_o, eby_o, delm_o, radius_o = @inferred uvbybeta(0.7901, 0.037, 0.039, 8)
    @test te_o ≈ 3985.833020767027
    @test mv_o ≈ 10.356215131980953
    @test eby_o ≈ 0.0027659846547315762
    @test delm_o ≈ 0.6944746189113468
    @test radius_o ≈ 0.2106831009861503
end

# Test airtovac is its inverse (it isn't true only around 2000, just avoid those values)
@testset "vactoair" begin
    @test @inferred(vactoair(2000)) ≈ 1999.3526230448367
    @test airtovac.(vactoair.(collect(1000:300:4000))) ≈ collect(1000:300:4000)
end

@testset "xyz" begin
    x, y, z, vx, vy, vz = @inferred(xyz([51200.5 + 64 / 86400], 2000))
    @test x  ≈ [0.5145687092402946]
    @test y  ≈ [-0.7696326261820777]
    @test z  ≈ [-0.33376880143026394]
    @test vx ≈ [0.014947267514081075]
    @test vy ≈ [0.008314838205475709]
    @test vz ≈ [0.003606857607574784]
end

@testset "ydn2md" begin
    @test ydn2md.(2016, [60, 234]) == [Date(2016, 02, 29), Date(2016, 08, 21)]
    @test @inferred(ymd2dn(ydn2md(2016, 60))) == 60
end

@testset "ymd2dn" begin
    @test ymd2dn.([Date(2015,3,5), Date(2016,3,5)]) == [64, 65]
    @test @inferred(ydn2md(2016, ymd2dn(Date(2016, 09, 16)))) == Date(2016, 09, 16)
end

@testset "zenpos" begin
    ra, dec = @inferred(zenpos(2.457514070138889e6, 45, 45))
    @test ra ≈ 1.9915758420649625
    @test dec ≈ 0.7853981633974483
    ra, dec = @inferred(zenpos(DateTime(2015, 11, 24, 13, 21), 43.16, -24.32, 4))
    @test ra ≈ 3.1232737646297757
    @test dec ≈ 0.7532841051607526
    ra, dec = @inferred(zenpos(jdcnv(2017, 01, 30, 04, 30), ten(35,0,42), ten(135,46,6)))
    @test ra ≈ 5.809762417608341
    @test dec ≈ 0.6110688599440813
end
