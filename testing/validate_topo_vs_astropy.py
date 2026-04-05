#!/usr/bin/env python3
"""
Validate the JS topocentric ICRF calculation against Astropy.

Takes geocentric EME2000 XYZ from the embedded OEM in artemis2_track.html,
computes topocentric RA/Dec two ways:
  1. Replica of the JS code (Lieske precession + IAU 1980 nutation + GAST)
  2. Astropy GCRS->ITRS machinery

Compares the results at a single observer location (Melbourne).
Expected result: <1 arcsec agreement (typically ~0.2" RMS).

Requires: numpy, astropy

Run from the testing/ directory:
    python3 validate_topo_vs_astropy.py
"""
import os
import numpy as np
import math
import re
from astropy.coordinates import EarthLocation, GCRS, ITRS, CartesianRepresentation
from astropy.time import Time
import astropy.units as u

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
HTML_PATH = os.path.join(SCRIPT_DIR, '..', 'artemis2_track.html')

R = math.pi / 180

# ── JS-equivalent functions ──

def gmst(jd):
    D = jd - 2451545.0
    T = D / 36525.0
    return (280.46061837 + 360.98564736629*D + 0.000387933*T*T - T*T*T/38710000.0) % 360

def nutation_vals(jd):
    T = (jd - 2451545.0) / 36525.0
    Om = (125.04 - 1934.136*T) * R
    Ls = (280.4665 + 36000.7698*T) * R
    Lm = (218.3165 + 481267.8813*T) * R
    dpsi = (-17.20*math.sin(Om) - 1.32*math.sin(2*Ls) - 0.23*math.sin(2*Lm) + 0.21*math.sin(2*Om)) / 3600 * R
    deps = (9.20*math.cos(Om) + 0.57*math.cos(2*Ls) + 0.10*math.cos(2*Lm) - 0.09*math.cos(2*Om)) / 3600 * R
    return dpsi, deps

def meanObliquity(jd):
    T = (jd - 2451545.0) / 36525.0
    return (23.439291111 - (46.8150 + 0.00059*T - 0.001813*T*T)*T/3600) * R

def gast_deg(jd):
    dpsi, _ = nutation_vals(jd)
    eps = meanObliquity(jd)
    eqEq = dpsi * math.cos(eps) * 180 / math.pi
    return (gmst(jd) + eqEq) % 360

def precessionMatrix(jd):
    T = (jd - 2451545.0) / 36525.0
    asR = math.pi / 648000.0
    zA = (2306.2181+1.39656*T-0.000139*T*T)*T + (0.30188-0.000344*T)*T*T + 0.017998*T*T*T
    zB = (2306.2181+1.39656*T-0.000139*T*T)*T + (1.09468+0.000066*T)*T*T + 0.018203*T*T*T
    th = (2004.3109-0.85330*T-0.000217*T*T)*T - (0.42665+0.000217*T)*T*T - 0.041833*T*T*T
    z1, z2, theta = zA*asR, zB*asR, th*asR
    cz1, sz1 = math.cos(z1), math.sin(z1)
    cz2, sz2 = math.cos(z2), math.sin(z2)
    ct, st = math.cos(theta), math.sin(theta)
    return np.array([
        [cz2*ct*cz1-sz2*sz1, -cz2*ct*sz1-sz2*cz1, -cz2*st],
        [sz2*ct*cz1+cz2*sz1, -sz2*ct*sz1+cz2*cz1, -sz2*st],
        [st*cz1, -st*sz1, ct]
    ])

def nutationMatrix(jd):
    dpsi, deps = nutation_vals(jd)
    eps = meanObliquity(jd)
    e1, e2 = eps, eps + deps
    ce1, se1 = math.cos(e1), math.sin(e1)
    ce2, se2 = math.cos(e2), math.sin(e2)
    cd, sd = math.cos(dpsi), math.sin(dpsi)
    return np.array([
        [cd, -sd*ce1, -sd*se1],
        [sd*ce2, cd*ce2*ce1+se2*se1, cd*ce2*se1-se2*ce1],
        [sd*se2, cd*se2*ce1-ce2*se1, cd*se2*se1+ce2*ce1]
    ])

def observerEME2000(lat, lon, alt_km, jd):
    a, f = 6378.137, 1/298.257223563
    e2 = 2*f - f*f
    phi, lam = lat*R, lon*R
    sp, cp = math.sin(phi), math.cos(phi)
    N = a / math.sqrt(1 - e2*sp*sp)
    xE = (N + alt_km)*cp*math.cos(lam)
    yE = (N + alt_km)*cp*math.sin(lam)
    zE = (N*(1-e2) + alt_km)*sp
    g = gast_deg(jd) * R
    cg, sg = math.cos(g), math.sin(g)
    xT = cg*xE - sg*yE
    yT = sg*xE + cg*yE
    zT = zE
    Nt = nutationMatrix(jd).T
    mean = Nt @ np.array([xT, yT, zT])
    Pt = precessionMatrix(jd).T
    return Pt @ mean

def topoRADec_js(x, y, z, jd, lat, lon, alt_km):
    obs = observerEME2000(lat, lon, alt_km, jd)
    dx, dy, dz = x - obs[0], y - obs[1], z - obs[2]
    r = math.sqrt(dx*dx + dy*dy + dz*dz)
    ra = math.atan2(dy, dx) / R % 360
    dec = math.asin(dz / r) / R
    return ra, dec, r

# ── Astropy calculation ──
def topoRADec_astropy(x_km, y_km, z_km, utc_str, lat, lon, alt_km):
    t = Time(utc_str, scale='utc')
    sc_gcrs = GCRS(
        CartesianRepresentation(x=x_km*u.km, y=y_km*u.km, z=z_km*u.km),
        obstime=t
    )
    loc = EarthLocation.from_geodetic(lon*u.deg, lat*u.deg, alt_km*1000*u.m)
    obs_itrs = loc.get_itrs(obstime=t)
    obs_gcrs = obs_itrs.transform_to(GCRS(obstime=t))
    dx = sc_gcrs.cartesian.x - obs_gcrs.cartesian.x
    dy = sc_gcrs.cartesian.y - obs_gcrs.cartesian.y
    dz = sc_gcrs.cartesian.z - obs_gcrs.cartesian.z
    r = np.sqrt(dx**2 + dy**2 + dz**2)
    ra = np.arctan2(dy, dx).to(u.deg).value % 360
    dec = np.arcsin(dz / r).to(u.deg).value
    return ra, dec, r.to(u.km).value


# ── Read OEM data from embedded HTML ──
oem_points = []
with open(HTML_PATH) as f:
    in_oem = False
    for line in f:
        if 'text/x-oem-data' in line:
            in_oem = True; continue
        if in_oem and '</script>' in line:
            break
        if not in_oem:
            continue
        m = re.match(r'(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)', line.strip())
        if m:
            oem_points.append((m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4))))

print(f"Loaded {len(oem_points)} OEM records from {HTML_PATH}")

# ── Test at various times ──
lat, lon, alt = -37.9039, 145.0224, 0.0
test_indices = [0, 50, 100, 200, 500, 800, 1000, 1500, 2000, 2500, 3000, len(oem_points)-1]

hdr = f"{'UTC':>25s}  {'JS_RA':>10s} {'AP_RA':>10s} {'dRA_as':>8s}  {'JS_Dec':>10s} {'AP_Dec':>10s} {'dDec_as':>8s}  {'sep_as':>7s}  {'range_km':>10s}"
print(hdr)
print('-' * len(hdr))

seps = []
for i in test_indices:
    if i >= len(oem_points):
        continue
    utc_str, x, y, z = oem_points[i]
    t = Time(utc_str, scale='utc')
    jd = t.jd

    ra_js, dec_js, rng_js = topoRADec_js(x, y, z, jd, lat, lon, alt)
    ra_ap, dec_ap, rng_ap = topoRADec_astropy(x, y, z, utc_str, lat, lon, alt)

    dra = (ra_js - ra_ap)
    if dra > 180: dra -= 360
    if dra < -180: dra += 360
    mean_dec = (dec_js + dec_ap) / 2
    dra_sky = dra * 3600 * math.cos(math.radians(mean_dec))
    ddec = (dec_js - dec_ap) * 3600
    sep = math.sqrt(dra_sky**2 + ddec**2)
    seps.append(sep)

    print(f"{utc_str:>25s}  {ra_js:10.5f} {ra_ap:10.5f} {dra_sky:+8.3f}  {dec_js:10.5f} {dec_ap:10.5f} {ddec:+8.3f}  {sep:7.3f}  {rng_js:10.1f}")

print(f"\nMax separation: {max(seps):.3f} arcsec")
print(f"RMS separation: {math.sqrt(sum(s**2 for s in seps)/len(seps)):.3f} arcsec")
