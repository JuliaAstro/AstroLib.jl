# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

"""
List of locations of North Magnetic Pole since 1590.

This is provided by World Magnetic Model
(https://www.ngdc.noaa.gov/geomag/data/poles/NP.xy).
"""
const POLELATLONG =
    try
        let
            local polelatlong, rows, floattype, temp
            polelatlong = readdlm(joinpath(dirname(@__FILE__),
                                           "..", "deps", "NP.xy"))
            rows = size(polelatlong, 1)
            floattype = typeof(polelatlong[1])
            temp = Dict{floattype, Tuple{floattype, floattype}}()
            for i = 1:rows
                merge!(temp, Dict(polelatlong[2rows + i]=>
                                  (polelatlong[rows + i],
                                   polelatlong[i])))
            end
            temp
        end
    catch
        warn("""Could not find file NP.xy, you will not be able to use
"geo2mag" and "mag2geo" functions.  Build again the package with
    Pkg.build("AstroLib")
then restart Julia session in order to fix this problem.""")
    end

const AU = 1.495978707e11 # Astronomical unit in meters
const J2000 = jdcnv(2000, 01, 01, 12) # 2451545
const JULIANYEAR = 365.25 # days in Julian year

# Constant used in ct2lst function, see Meeus, p.84.
const ct2lst_c  = (280.46061837, 360.98564736629, 0.000387933, 38710000.0)

# Used in "bprecess" and "jprecess".
const A_precess  = 1e-6*[-1.62557, -0.31919, -0.13843] # In radians
const A_dot_precess = 0.001*[1.244 , -1.579, -0.660] # In arc seconds per century

"""
List of observing sites.  The observatories have `Observatory` type.
"""
const observatories =
    Dict("mgio"=>Observatory("Mount Graham International Observatory",
                             (32,42,04.69), (-109,53,31.25), 3191.0, -7),
         "whitin"=>Observatory("Whitin Observatory, Wellesley College",
                               42.295, -71.305833, 32, -5),
         "fmo"=>Observatory("Fan Mountain Observatory",
                            (37,52,41), (-78,41,34), 556, -5),
         "lmo"=>Observatory("Leander McCormick Observatory",
                            (38,02,00), (-78,31,24), 264, -5),
         "holi"=>Observatory("Observatorium Hoher List (Universitaet Bonn) - Germany",
                             50.16276, 6.85, 541, 1),
         "ca"=>Observatory("Calar Alto Observatory",
                           (37,13,25), (-2,32,46.5), 2168, 1),
         "bgsuo"=>Observatory("Bowling Green State Univ Observatory",
                              (41,22,42), (-83,39,33), 225., -5),
         "irtf"=>Observatory("NASA Infrared Telescope Facility",
                             19.826218, -155.471999, 4168, -10),
         "rozhen"=>Observatory("National Astronomical Observatory Rozhen - Bulgaria",
                               (41,41,35), (24,44,38), 1759, 2),
         "bosque"=>Observatory("Estacion Astrofisica Bosque Alegre, Cordoba",
                               (-31,35,54), (-64,32,45), 1250, -3),
         "casleo"=>Observatory("Complejo Astronomico El Leoncito, San Juan",
                               (-31,47,57), (-69,18,00), 2552, -3),
         "saao"=>Observatory("South African Astronomical Observatory",
                             (-32,22,46), (20,48,38.5), 1798., 2),
         "lna"=>Observatory("Laboratorio Nacional de Astrofisica - Brazil",
                            (-22,32,04), -45.5825, 1864., -3),
         "oro"=>Observatory("Oak Ridge Observatory",
                            (42,30,18.94), (-71,33,29.32), 184., -5),
         "flwo"=>Observatory("Whipple Observatory",
                             (31,40,51.4), (-110,52,39), 2320., -7),
         "vbo"=>Observatory("Vainu Bappu Observatory",
                            12.57666, 78.8266, 725., 5.5),
         "lowell"=>Observatory("Lowell Observatory",
                               (35,05.8), (-111,32.1), 2198., -7),
         "apo"=>Observatory("Apache Point Observatory",
                            (32,46.8), (-105,49.2), 2798., -7),
         "loiano"=>Observatory("Bologna Astronomical Observatory, Loiano - Italy",
                               (44,15,33), (11,20,2), 785., 1),
         "ekar"=>Observatory("Mt. Ekar 182 cm. Telescope",
                             (45,50,54.92), (11,34,52.08), 1413.69, 1),
         "keck"=>Observatory("W. M. Keck Observatory",
                             (19,49.7), (-155,28.7), 4160., -10),
         "bao"=>Observatory("Beijing XingLong Observatory",
                            (40,23.6), (117,34.5), 950., 8),
         "bmo"=>Observatory("Black Moshannon Observatory",
                            (40,55.3), (-78,00.3), 738., -5),
         "nov"=>Observatory("National Observatory of Venezuela",
                            (8,47.4), (-70,52.0), 3610, -4),
         "mdm"=>Observatory("Michigan-Dartmouth-MIT Observatory",
                            (31,57.0), (-111,37.0), 1938.5, -7),
         "palomar"=>Observatory("The Hale Telescope",
                                (33,21,21.6), (-116,51,46.80), 1706., -8),
         # https://en.wikipedia.org/wiki/Tonantzintla_Observatory
         "tona"=>Observatory("Observatorio Astronomico Nacional, Tonantzintla",
                             (19,01,58), (-98,18,50), 2166., -6),
         "spm"=>Observatory("Observatorio Astronomico Nacional, San Pedro Martir",
                            (31,01,45), (-115,29,13), 2830., -7),
         "dao"=>Observatory("Dominion Astrophysical Observatory",
                            (48,31.3), (-123,25.0), 229., -8),
         "mtbigelow"=>Observatory("Catalina Observatory: 61 inch telescope",
                                  (32,25.0), (-110,43.9), 2510., -7),
         "lco"=>Observatory("Las Campanas Observatory",
                            (-29,0.2), (-70,42.1), 2282, -4),
         "mcdonald"=>Observatory("McDonald Observatory",
                                 30.6716667, -104.0216667, 2075, -6),
         "aao"=>Observatory("Anglo-Australian Observatory",
                            (-31,16,37.34), (149,3,57.91), 1164, 10),
         "sso"=>Observatory("Siding Spring Observatory",
                            (-31,16,24.10), (149,3,40.3), 1149, 10),
         "mso"=>Observatory("Mt. Stromlo Observatory",
                            (-35,19,14.34), (149,1,27.6), 767, 10),
         "lapalma"=>Observatory("Roque de los Muchachos, La Palma",
                                (28,45.5), (-17,52.8), 2327, 0),
         "cfht"=>Observatory("Canada-France-Hawaii Telescope",
                             (19,49.6), (-155,28.3), 4215., -10),
         "mmto"=>Observatory("MMT Observatory",
                             (31,41.3), (-110,53.1), 2600., -7),
         "lick"=>Observatory("Lick Observatory",
                             (37,20.6), (-121,38.2), 1290., -8),
         "eso"=>Observatory("European Southern Observatory",
                            (-29,15.4), (-70,43.8), 2347., -4),
         "ctio"=>Observatory("Cerro Tololo Interamerican Observatory",
                             -30.16527778, -70.815, 2215., -4),
         "kpno"=>Observatory("Kitt Peak National Observatory",
                             (31,57.8), (-111,36.0), 2120., -7))

"""
List of planets of the Solar System, from Mercury to Pluto.  The elements of the
list have `Planet` type.

Reference for most quantities is the Planetary Fact Sheet:
http://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
"""
const planets =
    Dict("mercury"=>Planet("mercury", 2439.7e3, 2439.7e3, 2439.7e3, 3.3011e23,
                           0.20563069, 57909050e3, 87.9691*86400),
         "venus"=>Planet("venus", 6051.8e3, 6051.8e3, 6051.8e3, 4.8675e24,
                         0.00677323, 108208000e3, 224.701*86400),
         "earth"=>Planet("earth", 6371e3, 6378.137e3,  6356.752e3, 5.97237e24,
                         0.0167086, 149598023e3, 365.256363004*86400),
         "mars"=>Planet("mars", 3389.5e3, 3396.2e3, 3376.2e3, 3.3011e23,
                        0.09341233, 57909050e3, 87.9691*86400),
         "jupiter"=>Planet("jupiter", 69911e3, 71492e3, 66854e3, 1898.19e24,
                           0.04839266, 5.20336301*AU, 4332.589*86400),
         "saturn"=>Planet("saturn", 58232e3, 60268e3, 54364e3, 568.34e23,
                          0.05415060, 9.53707032*AU, 10759.22*86400),
         "uranus"=>Planet("uranus", 25362e3, 25559e3, 24973e3, 86.813e24,
                          0.04716771, 19.19126393*AU, 30685.4*86400),
         "neptune"=>Planet("neptune", 24622e3, 24764e3, 24341e3, 102.413e24,
                           0.00858587, 30.06896348*AU, 60189*86400),
         "pluto"=>Planet("pluto", 1187e3, 1187e3, 1187e3, 0.01303e24,
                         0.24880766, 39.48168677*AU, 90560*86400))
