# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

# Test adstring
@test adstring((30.4, -1.23), truncate=true) ==
    adstring([30.4, -1.23], truncate=true) == " 02 01 35.9  -01 13 48"
@test adstring(19.19321, truncate=true) == "+19 11 35.5"
@test adstring(ten(36,24,15.015), -ten(8,24,36.0428), precision=3) ==
    " 02 25 37.0010  -08 24 36.043"
@test adstring([30.4, -15.63], [-1.23, 48.41], precision=2) ==
    [" 02 01 36.000  -01 13 48.00", " 22 57 28.800  +48 24 36.00"]
@test adstring([(-58, 724)]) == [" 20 08 00.0  +724 00 00"]

# Test airtovac
@test airtovac.([1234 6056.125]) ≈ [1234.0 6057.801930991426]
@test airtovac(2100) ≈ 2100.666421596007

# Test aitoff
@test aitoff([227.23, 130], [-8.890, -35]) ==
    ([-137.92196683723276, 115.17541338020645], [-11.772527357473054, -44.491889962090085])
@test aitoff([375], [2.437]) ==
    ([16.63760711611838],[2.712427279646118])
@test aitoff((227.23, -8.890)) ==
    (-137.92196683723276,-11.772527357473054)

# Test altaz2hadec
@test altaz2hadec(59.086111, 133.30806, 43.07833) ==
    altaz2hadec((59.086111, 133.30806), 43.07833) ==
    (336.68286017949157, 19.182449588316555)
@test altaz2hadec([15, 25, 35],
                  [25.12, 45.32, -20.3],
                  [-23.44, 45.0, 52.5]) ==
                      ([324.9881067314537, 256.7468302330436, 132.4919217875949],
                       [44.38225395397647, 48.542947077386664, 67.33061196497327])

# Test bprecess
let
    local ra, dec
    ra, dec = bprecess([ten(13, 42, 12.74)*15], [ten(8, 23, 17.69)],
                       reshape(100*[-15*0.0257, -0.090], 2, 1))
    @test ra  ≈ [204.93552515632123]
    @test dec ≈ [8.641287183886163]
    ra, dec = bprecess(82, 19)
    @test ra  ≈ 81.26467916346334
    @test dec ≈ 18.959495700195394
    ra, dec = bprecess([57], [23], 2024)
    @test ra  ≈ [56.26105898810067]
    @test dec ≈ [22.84693298145991]
    ra, dec = bprecess([57], [23], reshape([9, 86], 2, 1), parallax=[1], radvel=[4])
    @test ra  ≈ [56.25988479854577]
    @test dec ≈ [22.83493370392355]
    ra, dec = bprecess((-57, -23), 2024)
    @test ra  ≈  302.2593299643789
    @test dec ≈ -23.150089972802036
    ra, dec = bprecess((-57, -23), [9, 86], parallax=1, radvel=4)
    @test ra  ≈  302.2580376402947
    @test dec ≈ -23.16208183899836
end

# Test calz_unred
@test calz_unred.(collect(900:1000:9900), ones(Float64, 10), -0.1) ≈
    [1.0, 0.43189326452379095, 0.5203675483533704, 0.594996469192435, 0.6569506252451913,
     0.7080829505773865, 0.7502392743978797, 0.7861262388745882, 0.8151258710444882,
     0.8390325371659836]

#test ct2lst
@test ct2lst.(-76.72, -4, [DateTime(2008, 7, 30, 15, 53)]) ≈ [11.356505172312609]
@test ct2lst.(9, [jdcnv(2015, 11, 24, 12, 21)]) ≈ [17.159574059885927]

# Test daycnv with Gregorian Calendar in force.
@test daycnv(2440000.0) == DateTime(1968, 05, 23, 12)
# Test daycnv with Julian Calendar in force (same result as IDL AstroLib's
# daycnv).
@test daycnv(2000000.0) == DateTime(763, 09, 18, 12)
@test daycnv(0.0) == DateTime(-4713, 11, 24, 12)

# Test deredd
let
    by0, m0, c0, ub0 = deredd([0.5, -0.5], [0.2, 0.5], [1, 1], [1, 1], [0.1, 0.3])
    @test by0 ≈ [-0.3,0.5]
    @test m0  ≈ [1.165,1.0]
    @test c0  ≈ [0.905,1.0]
    @test ub0 ≈ [-0.665,0.3]
end


# Test eci2geo
let
    local lat, long, alt, jd
    lat, long, alt = eci2geo([0], [0], [0], [2452343])
    @test lat  ≈ [0]
    @test long ≈ [12.992783145436988]
    @test alt  ≈ [-6378.137]
    lat, long, alt = eci2geo((6978.137, 0, 0), jdcnv("2015-06-30T14:03:12.857"))
    @test lat  ≈ 0
    @test long ≈ 230.87301833205856
    @test alt  ≈ 600
    # Test `eci2geo' is the inverse of `geo2eci'
    jd = get_juldate()
    lat, long, alt = eci2geo(geo2eci(10, 10, 10, jd), jd)
    @test lat  ≈ 10
    @test long ≈ 10
    @test alt  ≈ 10
end

# Test eqpole
let
    local x, y
    x, y = eqpole([100], [35], southpole=true)
    @test x ≈ [-111.18287262822456]
    @test y ≈ [ -19.604540237028665]
    x, y = eqpole([80], [19])
    @test x ≈ [72.78853915267848]
    @test y ≈ [12.83458333897169]
end

# Test flux2mag
@test flux2mag.([1.5e-12, 8.7e-15, 4.4e-10]) ≈
    [8.459771852360795, 14.051201868453454, 2.291368308784527]
@test flux2mag(1) ≈ -21.1
@test flux2mag(5.2e-15) ≈ 14.609991640913002
@test flux2mag(5.2e-15, 15) ≈ 20.709991640913003
@test flux2mag(5.2e-15, ABwave=15) ≈ 27.423535345634598

# Test gal_uvw
let
    local u, v, w
    u, v, w = gal_uvw([ten(1,9,42.3)*15], [ten(61,32,49.5)], [627.89], [77.84],
                      [-321.4], [1e3/129], lsr=true)
    @test u ≈ [118.2110474553902]
    @test v ≈ [-466.4828898385057]
    @test w ≈ [88.16573278565097]
    u, v, w = gal_uvw(1, 2, 3, 4, 5, 6)
    @test u ≈  4.0228405867158745
    @test v ≈  3.7912174342038227
    @test w ≈ -3.1700191400725464
end

# Test geo2eci
let
    local x, y, z
    x, y, z = geo2eci([0], [0], [0], [2452343])
    @test x ≈ [6214.846433007192]
    @test y ≈ [-1433.9858454345972]
    @test z ≈ [0.0]
    x, y, z = geo2eci((0,0,0), jdcnv("2015-06-30T14:03:12.857"))
    @test x ≈ -4024.8671780315185
    @test y ≈ 4947.835465127513
    @test z ≈ 0.0
end

# Test geo2geodetic
let
    local lat, long, alt
    lat, long, alt = geo2geodetic([90], [0], [0], "Jupiter")
    @test lat  ≈   [90]
    @test long ≈    [0]
    @test alt  ≈ [4638]
    lat, long, alt = geo2geodetic((90, 0, 0))
    @test lat  ≈ 90
    @test long ≈  0
    @test alt  ≈ 21.38499999999931
    lat, long, alt = geo2geodetic((43.16, -24.32, 3.87), 8724.32, 8619.19)
    @test lat  ≈  43.849399515234516
    @test long ≈ -24.32
    @test alt  ≈  53.53354478670836
    lat, long, alt = geo2geodetic([43.16], [-24.32], [3.87], 8724.32, 8619.19)
    @test lat  ≈ [ 43.849399515234516]
    @test long ≈ [-24.32]
    @test alt  ≈ [ 53.53354478670836]
end

# Test geo2mag
let
    local lat, long
    lat, long = geo2mag(ten(35,0,42), ten(135,46,6), 2016)
    @test lat  ≈  36.86579228937769
    @test long ≈ -60.184060536651614
    lat, long = geo2mag([15], [24], 2016)
    @test lat  ≈ [ 11.452100529696096]
    @test long ≈ [-169.86030510727102]
end

# Test geodetic2geo
let
    local lat, long, alt
    lat, long, alt = geodetic2geo([90], [0], [0], "Jupiter")
    @test lat  ≈ [90]
    @test long == [0]
    @test alt  ≈ [-4638]
    lat, long, alt = geodetic2geo((90, 0, 0))
    @test lat ≈ 90
    @test long == 0
    @test alt ≈ -21.38499999999931
    lat, long, alt = geodetic2geo((43.16, -24.32, 3.87), 8724.32, 8619.19)
    @test lat ≈ 42.46772711708433
    @test long == -24.32
    @test alt ≈ -44.52902080669082
    lat, long, alt = geodetic2geo([43.16], [-24.32], [3.87], 8724.32, 8619.19)
    @test lat ≈ [42.46772711708433]
    @test long == [-24.32]
    @test alt ≈ [-44.52902080669082]
    # Test geodetic2geo is the inverse of geo2geodetic, within a certain
    # tolerance.
    lat, long, alt = geodetic2geo(geo2geodetic(67.2,13.4,1.2))
    @test lat ≈ 67.2 atol = 1e-8
    @test long == 13.4
    @test alt ≈ 1.2 atol = 1e-9
end

# Test get_date with mixed keywords.
@test get_date(DateTime(2001,09,25,14,56,14), old=true,timetag=true) ==
    get_date(2001,09,25,14,56,14, old=true,timetag=true) ==
    get_date("2001-09-25T14:56:14", old=true,timetag=true) ==
    "25/09/2001T14:56:14"
@test get_date(DateTime(2001,09,25,14,56,14)) ==
    get_date(2001,09,25,14,56,14) == get_date("2001-09-25T14:56:14") ==
    "2001-09-25"
@test get_date.([DateTime(2024), Date(2016, 3, 14)]) ==
    get_date.([Date(2024), "2016-03-14"]) ==
    get_date.(["2024-01", DateTime(2016, 3, 14)]) == ["2024-01-01", "2016-03-14"]

# Test gcirc.
@test gcirc.(0, [0,1,2], [1,2,3], [2,3,4], [3,4,5]) ≈
    [1.222450611061632, 2.500353926443337, 1.5892569925227757]
@test gcirc(0,  120, -43,   175, +22)     ≈  1.590442261600714
@test gcirc(1, (120, -43),  175, +22)     ≈  415908.56615322345
@test gcirc(2,  120, -43,  (175, +22))    ≈  296389.3666794745
@test gcirc(0, (120, -43), (175, +22))    ≈  1.590442261600714
@test gcirc.(1, [120], [-43],  175, +22)  ≈ [415908.56615322345]
@test gcirc.(2,  120, -43,  [175], [+22]) ≈ [296389.3666794745]
@test_throws ErrorException gcirc(3, 0, 0, 0, 0)

# Test hadec2altaz
let
    local alt1, az1, alt2, az2
    alt1, az1 = hadec2altaz([0], [11.978165], [ten(43,4,42)])
    @test alt1 ≈ [58.89983166666667]
    @test az1  ≈ [180.0]
    @test hadec2altaz((0, 11.978165), ten(43,4,42), ws=true)[2] ≈ 0.0
    alt1, az1 = 50, 20
    alt2, az2 = hadec2altaz(altaz2hadec(alt1, az1, 40), 40)
    @test alt1 ≈ alt2
    @test az1  ≈ az2
end

# Test helio_jd
@test helio_jd.([juldate(2016, 6, 15, 11, 40)], ten(20, 9, 7.8)*15, ten(37, 9, 7)) ≈
    [57554.98808289718]
@test helio_jd(1000, 23, 67, B1950=true) ≈ 999.9997659545342
@test helio_jd(2000, 12, 88, diff=true)  ≈ -167.24845957792076

# Test helio_rv
@test helio_rv(helio_jd(juldate(94, 10, 25, 17, 30),
                        ten(04, 38, 16)*15, ten(20, 41, 05)),
               46487.5303, 2.0563056, -6, 59.3) ≈ -62.965570109145034
@test helio_rv.([0.1, 0.9], 0, 1, 0, 100, 0.6, 45) ≈
    [-45.64994926111004, 89.7820347358485]

# Test jdcnv.
@test jdcnv(-4713, 11, 24, 12) ≈ 0.0
@test jdcnv(763, 09, 18, 12) == jdcnv("763-09-18T12") == 2000000.0
@test (jd=1234567.89; jdcnv(daycnv(jd)) == jd)
@test jdcnv.([DateTime(2016, 07, 31), "1969-07-20"]) ==
    jdcnv.([Date(2016, 07, 31), DateTime(1969, 07, 20)]) ==
    jdcnv.(["2016-07-31", Date(1969, 07, 20)])

# Test jprecess
let
    local ra, dec
    ra, dec = jprecess([ten(13, 39, 44.526)*15], [ten(8, 38, 28.63)],
                       reshape(100*[-15*0.0259, -0.093], 2, 1))
    @test ra  ≈ [205.5530845731372]
    @test dec ≈ [8.388247441904628]
    ra, dec = jprecess(82, 19)
    @test ra  ≈ 82.73568745151148
    @test dec ≈ 19.036972917272056
    ra, dec = jprecess([57], [23], 2024)
    @test ra  ≈ [57.74049975335702]
    @test dec ≈ [23.150053754297726]
    ra, dec = jprecess([57], [23], reshape([9, 86], 2, 1), parallax=[1], radvel=[4])
    @test ra  ≈ [57.74180294549785]
    @test dec ≈ [23.16200582079095]
    ra, dec = jprecess((-57, -23), 2024)
    @test ra  ≈ 303.73910971499015
    @test dec ≈ -22.846895476784482
    ra, dec = jprecess((-57, -23), [9, 86], parallax=1, radvel=4)
    @test ra  ≈ 303.7402950607101
    @test dec ≈ -22.834931625610313
end

# Test juldate with Gregorian Calendar in force.  This also makes sure precision
# of the result is high enough.  Note that "juldate(dt::DateTime) =
# Dates.datetime2julian(dt)-2.4e6" would not be precise.
@test juldate.([DateTime(2016, 1, 1, 8)]) ≈ [(57388.5 + 1//3)]

# Test juldate with Julian Calendar in force, for different centuries.  This
# also makes sure precision of the result is high enough.
@test juldate(1582, 10, 1, 20)    ≈ (-100843 + 1//3)
@test juldate.(["1000-01-01T20"]) ≈ [(-313692 + 1//3)]
@test juldate("100-10-25T20")     ≈ (-642119 + 1//3)
@test juldate(-4713, 1, 1, 12)    ≈ -2.4e6
@test juldate(2016, 06, 30, 00, 05, 53, 120) ≈
    jdcnv(2016, 06, 30, 00, 05, 53, 120) - 2.4e6
# Test daycnv and juldate together, with Gregorian Calendar in force.  Note that
# they are not expected to be one the inverse of the other during Julian
# Calendar.
@test (dt=DateTime(2016, 1, 1, 20, 45, 33, 457);
       daycnv(juldate(dt) + 2.4e6) == dt)

# Test kepler_solver
@test kepler_solver(8pi/3, 0.7)              ≈ 2.5085279492864223
@test kepler_solver.([pi/4, pi/6, 8pi/3], 0) ≈ [pi/4, pi/6, 2pi/3]
@test kepler_solver(3pi/2, 0.8)              ≈ -2.2119306096084457
@test kepler_solver(0, 1)                    ≈ 0.0
@test_throws ErrorException kepler_solver(pi, -0.5)
@test_throws ErrorException kepler_solver(pi,  1.5)

# Test lsf_rotate
let
    local vel, lsf
    vel, lsf = lsf_rotate(3, 90)
    @test length(vel) == length(lsf) == 61
    vel, lsf = lsf_rotate(5, 10)
    @test vel ≈ collect(-10.0:5.0:10.0)
    @test lsf ≈ [0.0, 0.556914447710896, 0.6933098861837907, 0.556914447710896, 0.0]
end

# Test mag2flux
@test mag2flux(4.83, 21.12)                         ≈ 4.1686938347033296e-11
@test mag2flux.([4.83], 21.12)                      ≈ [4.1686938347033296e-11]
@test flux2mag(mag2flux(15, ABwave=12.), ABwave=12) ≈ 15.0
@test mag2flux(8.3)                                 ≈ 1.7378008287493692e-12
@test mag2flux(8.3, 12)                             ≈ 7.58577575029182e-9
@test mag2flux(8.3, ABwave=12)                      ≈ 3.6244115683017193e-7

# Test mag2geo
let
    local lat, long
    lat, long = mag2geo(90, 0, 2016)
    @test lat  ≈ 86.395
    @test long ≈ -166.29000000000002
    lat, long = mag2geo([15], [24], 2016)
    @test lat  ≈ [11.702066965890157]
    @test long ≈ [-142.6357492442842]
    # Test geo2mag is approximately the inverse of mag2geo
    lat, long = geo2mag(mag2geo(12.34, 56.78, 2016)..., 2016)
    @test lat  ≈ 12.34
    @test long ≈ 56.78
end

# Test month_cnv
@test month_cnv.([" januavv  ", "SEPPES ", " aUgUsT", "la"]) == [1, 9, 8, -1]
@test month_cnv.([2, 12, 6], short=true, low=true) == ["feb", "dec", "jun"]
@test month_cnv(5, up=true) == "MAY"
@test (list=[1, 2, 3]; month_cnv.(month_cnv.(list)) == list)
@test (list=["July", "March", "November"]; month_cnv.(month_cnv.(list)) == list)

# Test moonpos
let
    local ra, dec, dis, lng, lat
    ra, dec, dis, lng, lat =
        moonpos(jdcnv(1992, 4, 12))
    @test ra  ≈ 134.68846854844108
    @test dec ≈ 13.768366630560255
    @test dis ≈ 368409.68481612665
    @test lng ≈ 133.16726428105378
    @test lat ≈ -3.2291264192144356
    ra, dec, dis, lng, lat =
        moonpos([2457521], radians=true)
    @test ra  ≈ [2.2587950290926178]
    @test dec ≈ [0.26183388011392217]
    @test dis ≈ [385634.68772395694]
    @test lng ≈ [2.232459255739553]
    @test lat ≈ [-0.059294466326164315]
end

# Test mphase
@test mphase.([2457520, 2457530, 2457650]) ≈
    [0.2781695910737857, 0.9969808583803166, 0.9580708477591693]

# Test nutate
let
    local long, obl
    long, obl = nutate(jdcnv(1987, 4, 10))
    @test long ≈ -3.787931077110755
    @test obl  ≈  9.442520698644401
    long, obl = nutate(2457521)
    @test long ≈ -4.401443629818089
    @test obl  ≈ -9.26823431959121
    long, obl = nutate([2457000, 2458000])
    @test long ≈ [ 4.327189321653877, -9.686089990639474]
    @test obl  ≈ [-9.507794266102866, -6.970768250588256]
end

# Test paczynski
@test paczynski(-1e-10) ≈  -1e10
@test paczynski(1e-1)   ≈  10.037461005722337
@test paczynski(-1)     ≈  -1.3416407864998738
@test paczynski(10)     ≈   1.0001922892047386
@test paczynski(-1e10)  ≈  -1

# Test planck_freq
@test planck_freq.([2000], [5000]) ≈ [6.1447146126144004e-30]

# Test planck_wave
@test planck_wave.([2000], [5000]) ≈ [8.127064833530511e-24]

# Test polrec
let
    local x, y
    x, y = polrec([1, 2, 3], [pi, pi/2.0, pi/4.0])
    @test x ≈ [-1.0, 0.0, 1.5*sqrt(2.0)]
    @test y ≈ [ 0.0, 2.0, 1.5*sqrt(2.0)]
    x, y = polrec((2, 135), degrees=true)
    @test x ≈ -sqrt(2)
    @test y ≈  sqrt(2)
end

# Test posang.
@test posang(1, ten(13, 25, 13.5), ten(54, 59, 17),
             ten(13, 23, 55.5), ten(54, 55, 31)) ≈ -108.46011246802047
@test posang.(0, [0,1,2], [1,2,3], [2,3,4], [3,4,5]) ≈
    [1.27896824717634, 1.6840484573313608, 0.2609280020139511]
@test posang(0,  120, -43,   175, +22)     ≈ -1.5842896165356724
@test posang(1, (120, -43),  175, +22)     ≈ 82.97831348792039
@test posang(2,  120, -43,  (175, +22))    ≈ 50.02816530382374
@test posang(0, (120, -43), (175, +22))    ≈ -1.5842896165356724
@test posang.(1, [120], [-43],  175, +22)  ≈ [82.97831348792039]
@test posang.(2,  120, -43,  [175], [+22]) ≈ [50.02816530382374]
@test_throws ErrorException posang(3, 0, 0, 0, 0)

# Test precess
let
    local ra1, dec1, ra2, dec2
    ra1, dec1 = precess((ten(2,31,46.3)*15, ten(89,15,50.6)), 2000, 1985)
    @test ra1  ≈ 34.09470328718033
    @test dec1 ≈ 89.19647174928589
    ra2, dec2 = precess([ten(21, 59, 33.053)*15], [ten(-56, 59, 33.053)], 1950, 1975, FK4=true)
    @test ra2  ≈ [330.3144305418865]
    @test dec2 ≈ [-56.87186126487889]
end

# Test precess_xyz
let
    local x1 ,y1, z1, x2, y2, z2
    x1, y1, z1 = precess_xyz((1.2, 2.3, 1.7), 2000, 2050)
    @test x1 ≈ 1.165933061423247
    @test y1 ≈ 2.313228746401996
    @test z1 ≈ 1.7057470102860104
    x2, y2, z2 = precess_xyz([0.7, -2.4], [3.3, 6.6], [0, 4], 2000, 2016)
    @test x2 ≈ [0.688187142071843,   -2.429815562246262]
    @test y2 ≈ [3.3024835038223532,   6.591359330834213]
    @test z2 ≈ [0.001079105285993004, 3.9962455511755794]
end

# Test premat
@test premat(1967, 1982, FK4=true) ≈
    [0.9999933170034135    -0.0033529069683496567 -0.0014573823699636742;
     0.00335290696825777    0.9999943789886484    -2.443304965138481e-6 ;
     0.0014573823701750721 -2.4431788671274868e-6  0.9999989380147651   ]
@test premat(1995, 2003) ≈
    [ 0.9999980977132219    -0.0017889257711428855 -0.0007773766929507687;
      0.0017889257711354528  0.9999983998707707    -6.953448226403318e-7 ;
      0.0007773766929678732 -6.953257000046125e-7   0.9999996978424512   ]

# Test radec
@test radec(15.90, -0.85) == (1.0, 3.0, 36.0, -0.0, 51.0, 0.0)
@test radec(-0.85,15.9) == (23.0,56.0,36.0,15.0,54.0,0.0)
@test radec(-20,4,hours=true) == (4.0,0.0,0.0,4.0,0.0,0.0)
@test radec([15.90, -0.85], [-0.85,15.9]) ==
    ([1.0, 23.0], [3.0, 56.0], [36.0, 36.0],
     [-0.0, 15.0], [51.0, 54.0], [0.0, 0.0])

# Test recpol
let
    r = a = zeros(Float64, 3)
    r, a = recpol([0, sqrt(2.0), 2.0*sqrt(3.0)],
                  [0, sqrt(2.0), 2.0])
    @test r ≈ [0.0,  2.0,  4.0]
    @test a ≈ [0.0, pi/4.0, pi/6.0]
    r, a = recpol(1, 1)
    @test r ≈ sqrt(2.0)
    @test a ≈ pi/4.0
    # Test polrec is the inverse of recpol
    xi, yi, x, y, = 6.3, -2.7, 0.0, 0.0
    x, y = polrec(recpol((xi, yi), degrees=true), degrees=true)
    @test x ≈ xi
    @test y ≈ yi
end

# Test rhotheta
let
    local ρ, θ
    ρ, θ = rhotheta(41.623, 1934.008, 0.2763, 0.907, 59.025, 23.717, 219.907, 1980)
    @test ρ ≈ 0.41101776646245836
    @test θ ≈ 318.4242564860495
end

# Test "sixty".  Test also it's the reverse of ten.
@test sixty(-51.36)                    ≈ [-51.0, 21.0, 36.0]
@test ten(sixty(-0.10934835545824395)) ≈ -0.10934835545824395
@test sixty(1)                         ≈ [1.0, 0.0, 0.0]

# Test sphdist.
@test sphdist.([0,1,2], [1,2,3], [2,3,4], [3,4,5]) ≈
    [1.222450611061632, 2.500353926443337, 1.5892569925227762]
@test sphdist(120, -43, 175, +22)      ≈  1.5904422616007134
@test sphdist.([120], [-43], 175, +22) ≈ [1.5904422616007134]
@test sphdist.(120, -43, [175], [+22]) ≈ [1.5904422616007134]

# Test sunpos
let
    local ra, dec, lon, obl
    ra, dec, lon, obl = sunpos(jdcnv(1982, 5, 1))
    @test ra  ≈ 37.88589057369026
    @test dec ≈ 14.909699471099517
    @test lon ≈ 2309.6312912217463
    @test obl ≈ 1343.0612563977593
    ra, dec, lon, obl = sunpos(jdcnv.([DateTime(2016, 5, 10)]), radians=true)
    @test ra  ≈ [0.8259691339090751]
    @test dec ≈ [0.3085047454107549]
    @test lon ≈ [49.77773359512005 ]
    @test obl ≈ [23.434647165246304]
    ra, dec, lon, obl = sunpos([2457531])
    @test ra  ≈ [59.71655864208797]
    @test dec ≈ [20.52127006818727]
    @test lon ≈ [3542.279299068626]
    @test obl ≈ [1342.7064622183311]
end

# Test "ten" and "tenv".  Always make sure string and numerical inputs are
# consistent (IDL implementation of "ten" is not).
@test ten(0, -23, 34) == ten((0, -23, 34)) == ten([0, -23, 34]) ==
    ten(" : 0 -23 :: 34") == -0.37388888888888894
@test ten(-0.0, 60) == ten((-0.0, 60)) == ten([-0.0, 60]) ==
    ten("-0.0 60") == -1.0
@test ten(-5, -60, -3600) == ten((-5, -60, -3600)) ==
    ten([-5, -60, -3600]) == ten("  -5: :-60: -3600") == -3.0
@test ten("") == 0.0
@test ten.([0, -0.0, -5], [-23, 60, -60], [34, 0, -3600]) ==
    ten.([(0, -23,34), ":-0.0:60", (-5, -60, -3600)]) ==
    ten.(["0   -23 :: 34", (-0.0, 60), " -5:-60: -3600"]) ==
    [-0.37388888888888894, -1.0, -3.0]
@test ten.([12.0, -0.0], [24, 30]) == ten.([" 12::24", " -0:30: "]) == [12.4, -0.5]

# Test kepler_solver
@test trueanom(8pi/3, 0.7)              ≈ 2.6657104039293764
@test trueanom.([pi/4, pi/6, 8pi/3], 0) ≈ [pi/4, pi/6, 2pi/3]
@test trueanom(3pi/2, 0.8)              ≈ -2.498091544796509
@test trueanom(0.1, 1)                  ≈ pi
@test_throws ErrorException trueanom(pi, -0.5)
@test_throws DomainError trueanom(pi,  1.5)

# Test vactoair and that airtovac is its inverse (it isn't true only around
# 2000, just avoid those values)
@test vactoair.([2000]) ≈ [1999.3526230448367]
@test airtovac.(vactoair.(collect(1000:300:4000))) ≈ collect(1000:300:4000)

# Test xyz
let
    local x, y, z, vx, vy, vz
    x, y, z, vx, vy, vz = xyz([51200.5 + 64./86400.], 2000)
    @test x  ≈ [0.5145687092402946]
    @test y  ≈ [-0.7696326261820777]
    @test z  ≈ [-0.33376880143026394]
    @test vx ≈ [0.014947267514081075]
    @test vy ≈ [0.008314838205475709]
    @test vz ≈ [0.003606857607574784]
end

# Test ydn2md.
@test ydn2md.(2016, [60, 234]) == [Date(2016, 02, 29), Date(2016, 08, 21)]
@test ymd2dn(ydn2md(2016, 60)) == 60

# Test ymd2dn
@test ymd2dn.([Date(2015,3,5), Date(2016,3,5)]) == [64, 65]
@test ydn2md(2016, ymd2dn(Date(2016, 09, 16))) == Date(2016, 09, 16)
