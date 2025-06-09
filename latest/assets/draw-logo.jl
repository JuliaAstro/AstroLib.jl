#!/usr/bin/env julia

using Luxor

function ellipse2c(focus1::Point, focus2::Point, k; stepvalue=pi/100)
    # k is sum of distances to focii of any points on the ellipse)
    a = k/2  # Major axis
    cpoint = midpoint(focus1, focus2)
    dc = norm(focus1, cpoint) # Distance of the foci to the center
    b = sqrt(abs(a^2 - dc^2)) # Minor axis
    phi = slope(O, focus2) # Angle between the major axis and the x-axis
    drawing = false
    for t in 0:stepvalue:2pi
        xt = cpoint.x + a * cos(t) * cos(phi) - b * sin(t) * sin(phi)
        yt = cpoint.y + a * cos(t) * sin(phi) + b * sin(t) * cos(phi)
        if drawing
            line(Point(xt, yt))
        else
            move(Point(xt, yt))
            drawing = true
        end
    end
    closepath()
    stroke()
end

"""
Draw a logo with ellipse and colored circles 
"""
function draw_logo(fname)
    Drawing(500, 500, fname)
    setline(3)
    origin()
    cpoints = ngon(O, 200, 3, vertices=true)

    ellipse2c(O, cpoints[1], 280)
    ellipse2c(O, cpoints[2], 280)
    ellipse2c(O, cpoints[3], 280)

    darker_purple = (0.584, 0.345, 0.698)
    lighter_purple  = (0.667, 0.475, 0.757)
    darker_green  = (0.22, 0.596, 0.149)
    lighter_green  = (0.376, 0.678, 0.318)
    darker_red  = (0.796, 0.235, 0.2)
    lighter_red  = (0.835, 0.388, 0.361)

    ballradius = 40
    sethue(darker_red...)
    circle(cpoints[1], 0.75ballradius, :fill)
    sethue(lighter_red...)
    circle(cpoints[1], 0.65ballradius, :fill)

    sethue(darker_green...)
    circle(cpoints[2], 0.75ballradius, :fill)
    sethue(lighter_green...)
    circle(cpoints[2], 0.65ballradius, :fill)

    sethue(darker_purple...)
    circle(cpoints[3], 0.75ballradius, :fill)
    sethue(lighter_purple...)
    circle(cpoints[3], 0.65ballradius, :fill)

    sethue("gold")
    star(O, 35, 12, 0.7, pi/3, :fill)
    sethue("yellow")
    star(O, 35, 15, 0.3, 0, :fill)
    finish()
    preview()
end

draw_logo("/tmp/logo.png")
