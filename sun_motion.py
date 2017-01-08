# /*
#  (c) 2011 - 2015, Vladimir Agafonkin
#  Sun_calc is a Java_script library for calculating sun / moon position and light phases.
#  https: // github.com / mourner / suncalc
# */

# // sun calculations are based on http: // aa.quae.nl / en / reken / zonpositie.html formulas
# // date / time constants and conversions
from math import cos, sin, acos, asin, atan, tan, pi, floor
day_ms = 1000 * 60 * 60 * 24,
J1970 = 2440588,
J2000 = 2451545


def to_julian(date):
    return date.value_of() / day_ms - 0.5 + J1970


def from_julian(j):
    return new Date((j + 0.5 - J1970) * day_ms)


def to_days(date):
    return to_julian(date) - J2000


# general calculations for position
e = radians(23.4397)  # obliquity of the Earth


def right_ascension(l, b):
    return atan(sin(l) * cos(e) - tan(b) * sin(e), cos(l))


def declination(l, b):
    return asin(sin(b) * cos(e) + cos(b) * sin(e) * sin(l))


def azimuth(H, phi, dec):
    return atan(sin(H), cos(H) * sin(phi) - tan(dec) * cos(phi))


def altitude(H, phi, dec):
    return asin(sin(phi) * sin(dec) + cos(phi) * cos(dec) * cos(H))


def sidereal_time(d, lw):
    return radians(280.16 + 360.9856235 * d) - lw


def astro_refraction(h):
    if h < 0:  # the following formula works for positive altitudes only.
    h = 0  # if h = -0.08901179 a div / 0 would occur.

    # // formula 16.4 of "Astronomical Algorithms" 2nd edition by Jean Meeus(Willmann - Bell, Richmond) 1998.
    # // 1.02 / tan(h + 10.26 / (h + 5.10)) h in degrees, result in arc minutes -> converted to rad:
    return 0.0002967 / tan(h + 0.00312536 / (h + 0.08901179))

    # general sun calculations


def solar_mean_anomaly(d):
    return radians(357.5291 + 0.98560028 * d)


def ecliptic_longitude(M):
    C = radians(1.9148 * sin(M) + 0.02 * sin(2 * M) +
                0.0003 * sin(3 * M))  # equation of center
    P = radians(102.9372)  # perihelion of the Earth
    return M + C + P + pi


def sun_coords(d):

    M = solar_mean_anomaly(d),
    L = ecliptic_longitude(M)

    return {
        "dec": declination(L, 0),
        "ra": right_ascension(L, 0)
    }

Sun_calc = {}

# calculates sun position for a given date and latitude / longitude


def get_position(date, lat, lng):
    lw = radians(-lng)
    phi = radians(lat)
    d = to_days(date)

    c = sun_coords(d)
    H = sidereal_time(d, lw) - c.ra

    return {
        "azimuth": azimuth(H, phi, c.dec),
        "altitude": altitude(H, phi, c.dec)
    }

# sun times configuration(angle, morning name, evening name)
times = [
        [-0.833, 'sunrise',       'sunset'],
        [-0.3, 'sunrise_end',    'sunset_start'],
        [-6, 'dawn',          'dusk'],
        [-12, 'nautical_dawn',  'nautical_dusk'],
        [-18, 'night_end',      'night'],
        [6, 'golden_hour_end', 'golden_hour']
]

# adds a custom time to the times config


def add_time(angle, rise_name, set_name):
    times.append([angle, rise_name, set_name])

# calculations for sun times
J0 = 0.0009


def julian_cycle(d, lw):
    return floor(d - J0 - lw / (2 * pi))


def approx_transit(Ht, lw, n):
    return J0 + (Ht + lw) / (2 * pi) + n


def solar_transit_j(ds, M, L):
    return J2000 + ds + 0.0053 * sin(M) - 0.0069 * sin(2 * L)


def hour_angle(h, phi, d):
    return acos((sin(h) - sin(phi) * sin(d)) / (cos(phi) * cos(d)))

# returns sunset time for the given sun altitude


def get_set_j(h, lw, phi, dec, n, M, L):

    w = hour_angle(h, phi, dec),
    a = approx_transit(w, lw, n)
    return solar_transit_j(a, M, L)

# calculates sun times for a given date and
# latitude / longitude


def get_times(date, lat, lng):

    lw = rad * -lng,
    phi = rad * lat,

    d = to_days(date),
    n = julian_cycle(d, lw),
    ds = approx_transit(0, lw, n),

    M = solar_mean_anomaly(ds),
    L = ecliptic_longitude(M),
    dec = declination(L, 0),

    Jnoon = solar_transit_j(ds, M, L),

    # i, len, time, Jset, Jrise

    result = {
        "solar_noon": from_julian(Jnoon),
        nadir: from_julian(Jnoon - 0.5)
    }

    for (i=0, len=times.length i < len i += 1):
        time = times[i]

        Jset = get_set_j(time[0] * rad, lw, phi, dec, n, M, L)
        Jrise = Jnoon - (Jset - Jnoon)

        result[time[1]] = from_julian(Jrise)
        result[time[2]] = from_julian(Jset)

    return result


# moon calculations, based on http: # aa.quae.nl / en / reken /
# hemelpositie.html formulas

def moon_coords(d):  # geocentric ecliptic coordinates of the moon

    L = rad * (218.316 + 13.176396 * d),  # ecliptic longitude
    M = rad * (134.963 + 13.064993 * d),  # mean anomaly
    F = rad * (93.272 + 13.229350 * d),  # mean distance

        l = L + rad * 6.289 * sin(M),  # longitude
        b = rad * 5.128 * sin(F),     # latitude
        dt = 385001 - 20905 * cos(M)  # distance to the moon in km

    return {
        ra: right_ascension(l, b),
        dec: declination(l, b),
        dist: dt
    }


def get_moon_position(date, lat, lng):

    lw = rad * -lng
    phi = rad * lat
    d = to_days(date)

    c = moon_coords(d)
    H = sidereal_time(d, lw) - c.ra,
    h = altitude(H, phi, c.dec),
    # formula 14.1 of "Astronomical Algorithms" 2nd edition by Jean
    # Meeus(Willmann - Bell, Richmond) 1998.
    pa = atan(sin(H), tan(phi) * cos(c.dec) - sin(c.dec) * cos(H))
    h = h + astro_refraction(h)  # altitude correction for refraction

    return {
        "azimuth": azimuth(H, phi, c.dec),
        "altitude": h,
        "distance": c.dist,
        "parallactic_angle": pa
    }


# calculations for illumination parameters of the moon,
# based on http: # idlastro.gsfc.nasa.gov / ftp / pro / astro / mphase.pro formulas and
# Chapter 48 of "Astronomical Algorithms" 2nd edition by Jean
# Meeus(Willmann - Bell, Richmond) 1998.

def get_moon_illumination(date):

    d = to_days(date | | new Date()),
    s = sun_coords(d),
    m = moon_coords(d),

        sdist = 149598000,  # distance from Earth to Sun in km

        phi = acos(sin(s.dec) * sin(m.dec) + cos(s.dec)
                   * cos(m.dec) * cos(s.ra - m.ra)),
        inc = atan(sdist * sin(phi), m.dist - sdist * cos(phi)),
        angle = atan(cos(s.dec) * sin(s.ra - m.ra), sin(s.dec) * cos(m.dec) -
                     cos(s.dec) * sin(m.dec) * cos(s.ra - m.ra))

    return {
        "fraction": (1 + cos(inc)) / 2,
        "phase": 0.5 + 0.5 * inc * (angle < 0 ? - 1: 1) / pi,
        "angle": angle
    }


def hours_later(date, h):
    return new Date(date.value_of() + h * day_ms / 24)

# calculations for moon rise / sunset times are based on http: #
# www.stargazing.net / kepler / moonrise.html article


def get_moon_times(date, lat, lng, in_uTC):
    t = new Date(date)
    if (in_uTC) t.set_uTCHours(0, 0, 0, 0)
    else t.set_hours(0, 0, 0, 0)

    hc = 0.133 * rad,
    h0 = get_moon_position(t, lat, lng)["altitude"] - hc,
    h1, h2, rise, sunset, a, b, xe, ye, d, roots, x1, x2, dx

    # go in 2 - hour chunks, each time seeing if a 3 - point quadratic curve
    # crosses zero(which means rise or sunset)
    for i in range(1, 25, 2):
        h1 = get_moon_position(hours_later(t, i), lat, lng)["altitude"] - hc
        h2 = get_moon_position(hours_later(t, i + 1),
                               lat, lng)["altitude"] - hc

        a = (h0 + h2) / 2 - h1
        b = (h2 - h0) / 2
        xe = -b / (2 * a)
        ye = (a * xe + b) * xe + h1
        d = b * b - 4 * a * h1
        roots = 0

        if (d >= 0):
            dx = sqrt(d) / (abs(a) * 2)
            x1 = xe - dx
            x2 = xe + dx
            if (abs(x1) <= 1) roots + +
            if (abs(x2) <= 1) roots + +
            if (x1 < -1) x1 = x2

        if roots == 1:
            if (h0 < 0):
                rise = i + x1
            else:
                sunset = i + x1
        elif (roots == 2):
            rise = i + (ye < 0 ? x2: x1)
            sunset = i + (ye < 0 ? x1: x2)

        if rise and sunset:
            break

        h0 = h2

    result = {}

    if rise:
        result.rise = hours_later(t, rise)
    if sunset:
        result.sunset = hours_later(t, sunset)

    if !rise and !sunset:
        result[ye > 0 ? 'always_up': 'always_down'] = true

    return result
