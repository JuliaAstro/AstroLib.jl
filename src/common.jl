# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

using DelimitedFiles

"""
List of locations of North Magnetic Pole since 1590.

This is provided by World Magnetic Model
(https://www.ngdc.noaa.gov/geomag/data/poles/NP.xy).
"""
const POLELATLONG =
    let
        local polelatlong, rows, floattype, temp
        polelatlong = readdlm(joinpath(@__DIR__,
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

const AU = 149_597_870_700 # Astronomical unit in meters
const J2000 = Int(jdcnv(2000, 01, 01, 12)) # 2451545
const JULIANYEAR = 365.25 # days in one Julian year
const JULIANCENTURY = 36_525 # days in one Julian century

# Constant used in ct2lst function, see Meeus, p.84.
const ct2lst_c  = (280.46061837, 360.98564736629, 0.000387933, 38710000.0)

# Used in "bprecess" and "jprecess".
const A_precess  = SVector(-1.62557, -0.31919, -0.13843) ./ 1000000 # In radians
const A_dot_precess = SVector(1.244 , -1.579, -0.660) ./ 1000 # In arc seconds per century

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
                             (31,57.8), (-111,36.0), 2120., -7),
         # https://en.wikipedia.org/wiki/Pine_Bluff_Observatory
         "pbo"=>Observatory("Pine Bluff Observatory",
                            43.0777, -89.6717, 362, -6))

"""
List of planets of the Solar System, from Mercury to Pluto.  The elements of the
list have `Planet` type.

Reference for most quantities is the Planetary Fact Sheet:
http://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
and the Keplerian Elements for Approximate Positions of the
Major Planets: https://ssd.jpl.nasa.gov/txt/p_elem_t1.txt
"""
const planets =
    Dict("mercury"=>Planet("mercury", 2439.7e3, 2439.7e3, 2439.7e3, 3.3011e23,
                           0.20563593, 0.38709927 *AU, 87.9691*86400,
                           7.00497902, 48.33076593, 77.45779628, 252.25032350),
         "venus"=>Planet("venus", 6051.8e3, 6051.8e3, 6051.8e3, 4.8675e24,
                         0.00677672, 0.72333566 * AU, 224.701*86400,
                         3.39467605, 76.67984255, 131.60246718, 181.97909950),
         "earth"=>Planet("earth", 6371e3, 6378.137e3,  6356.752e3, 5.97237e24,
                         0.01671123, 1.00000261 * AU, 365.256363004*86400,
                         -0.00001531, 0, 102.93768193, 100.46457166),
         "mars"=>Planet("mars", 3389.5e3, 3396.2e3, 3376.2e3, 3.3011e23,
                        0.09339410, 1.52371034 * AU, 87.9691*86400,
                        1.84969142, 49.55953891, -23.94362959, -4.55343205),
         "jupiter"=>Planet("jupiter", 69911e3, 71492e3, 66854e3, 1898.19e24,
                           0.04838624, 5.20288700 * AU, 4332.589*86400,
                           1.30439695, 100.47390909, 14.72847983, 34.39644051),
         "saturn"=>Planet("saturn", 58232e3, 60268e3, 54364e3, 568.34e24,
                          0.05386179, 9.53667594 * AU, 10759.22*86400,
                          2.48599187, 113.66242448, 92.59887831, 49.95424423),
         "uranus"=>Planet("uranus", 25362e3, 25559e3, 24973e3, 86.813e24,
                          0.04725744, 19.18916464 * AU, 30685.4*86400,
                          0.77263783, 74.01692503, 170.95427630, 313.23810451),
         "neptune"=>Planet("neptune", 24622e3, 24764e3, 24341e3, 102.413e24,
                           0.00859048, 30.06992276 * AU, 60189.0*86400,
                           1.77004347, 131.78422574, 44.96476227, -55.12002969),
         "pluto"=>Planet("pluto", 1187e3, 1187e3, 1187e3, 0.01303e24,
                         0.24882730, 39.48211675 * AU, 90560.0*86400,
                         17.14001206, 110.30393684, 224.06891629, 238.92903833))

if !isdefined(Base, :sincos)
    sincos(x::Real) = (sin(x), cos(x))
end
