#!/usr/bin/env python3
"""
Detect scheduled burns in Artemis II trajectory data.

Method:
  1. Convert geocentric (RA, Dec, Range) to XYZ position vectors
  2. Compute velocity via Savitzky-Golay differentiation (smooth, noise-resistant)
  3. Compute acceleration the same way
  4. Subtract expected gravitational acceleration (Earth + Moon two-body)
  5. The residual acceleration = non-gravitational force = BURNS
  6. Spikes in |residual| flag burn times

Input: RA/Dec text file exported from the Artemis II web tracker.
"""

import sys
from datetime import datetime, timedelta
import math
import numpy as np

# ---------- parse ----------
def parse_file(path):
    records = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(('Artemis', 'Generated', 'AEST', '─')):
                continue
            parts = line.split()
            if len(parts) < 16:
                continue
            try:
                utc = datetime.fromisoformat(parts[3])
            except (ValueError, IndexError):
                continue
            try:
                records.append({
                    'utc': utc,
                    'ra': float(parts[4]),
                    'dec': float(parts[5]),
                    'dra': float(parts[12]),
                    'ddec': float(parts[13]),
                    'range': float(parts[14]),
                    'mot': float(parts[15]),
                })
            except (ValueError, IndexError):
                continue
    return records

# ---------- Moon position (simplified Meeus) ----------
def moon_geocentric_xyz(jd):
    """Approximate Moon geocentric XYZ (EME2000, km). Good to ~0.5°/1000 km."""
    T = (jd - 2451545.0) / 36525.0
    # Mean elements (degrees)
    Lp = 218.3165 + 481267.8813 * T  # mean longitude
    D  = 297.8502 + 445267.1115 * T  # mean elongation
    M  = 357.5291 +  35999.0503 * T  # Sun mean anomaly
    Mp = 134.9634 + 477198.8676 * T  # Moon mean anomaly
    F  = 93.2720  + 483202.0175 * T  # argument of latitude

    Lp, D, M, Mp, F = [x % 360 for x in [Lp, D, M, Mp, F]]
    dr = math.pi / 180

    # Longitude and latitude perturbations (largest terms)
    lon = Lp + 6.289 * math.sin(Mp*dr) \
            + 1.274 * math.sin((2*D - Mp)*dr) \
            + 0.658 * math.sin(2*D*dr) \
            + 0.214 * math.sin(2*Mp*dr) \
            - 0.186 * math.sin(M*dr) \
            - 0.114 * math.sin(2*F*dr)
    lat = 5.128 * math.sin(F*dr) \
        + 0.281 * math.sin((Mp + F)*dr) \
        + 0.278 * math.sin((Mp - F)*dr)

    # Distance (km)
    dist = 385001 - 20905 * math.cos(Mp*dr) \
                  -  3699 * math.cos((2*D - Mp)*dr) \
                  -  2956 * math.cos(2*D*dr)

    # Ecliptic to equatorial (obliquity ≈ 23.4393°)
    eps = 23.4393 * dr
    lon_r = lon * dr
    lat_r = lat * dr

    # Ecliptic XYZ
    xe = dist * math.cos(lat_r) * math.cos(lon_r)
    ye = dist * math.cos(lat_r) * math.sin(lon_r)
    ze = dist * math.sin(lat_r)

    # Rotate to equatorial
    x = xe
    y = ye * math.cos(eps) - ze * math.sin(eps)
    z = ye * math.sin(eps) + ze * math.cos(eps)
    return np.array([x, y, z])

# ---------- Savitzky-Golay differentiation ----------
def savgol_coeffs(window_length, polyorder, deriv=0):
    """Compute Savitzky-Golay filter coefficients for given derivative."""
    half = window_length // 2
    # Build Vandermonde matrix
    x = np.arange(-half, half + 1, dtype=float)
    A = np.zeros((window_length, polyorder + 1))
    for j in range(polyorder + 1):
        A[:, j] = x ** j
    # Least squares: coeffs = (A^T A)^-1 A^T, row for deriv
    ATA = A.T @ A
    ATAinv = np.linalg.inv(ATA)
    coeffs = (ATAinv @ A.T)[deriv] * math.factorial(deriv)
    return coeffs

def savgol_filter(y, window_length, polyorder, deriv=0, dt=1.0):
    """Apply Savitzky-Golay filter/differentiation."""
    coeffs = savgol_coeffs(window_length, polyorder, deriv)
    half = window_length // 2
    n = len(y)
    out = np.zeros(n)
    # Convolve (with edge padding by reflection)
    padded = np.concatenate([y[half:0:-1], y, y[-2:-half-2:-1]])
    for i in range(n):
        out[i] = np.dot(coeffs, padded[i:i+window_length])
    return out / (dt ** deriv)

# ---------- main analysis ----------
def main():
    path = sys.argv[1] if len(sys.argv) > 1 else '/Users/karl/Downloads/artemis2_radec-2.txt'
    print(f"Reading: {path}")
    print(f"(Input: RA/Dec text export from Artemis II web tracker)\n")
    records = parse_file(path)
    n = len(records)
    print(f"Parsed {n} records at 10-min intervals")
    print(f"Time span: {records[0]['utc']} → {records[-1]['utc']}")

    dt = 600.0  # seconds

    # Convert to XYZ (geocentric EME2000)
    xyz = np.zeros((n, 3))
    for i, r in enumerate(records):
        ra_r = math.radians(r['ra'])
        dec_r = math.radians(r['dec'])
        d = r['range']
        xyz[i] = [d * math.cos(dec_r) * math.cos(ra_r),
                  d * math.cos(dec_r) * math.sin(ra_r),
                  d * math.sin(dec_r)]

    # Compute velocity and acceleration via S-G filter
    # Window = 11 samples (110 min), polynomial order 5
    # This is wide enough to smooth noise but narrow enough to resolve burns
    win = 11
    poly = 5
    vel = np.zeros((n, 3))
    acc_obs = np.zeros((n, 3))
    for axis in range(3):
        vel[:, axis] = savgol_filter(xyz[:, axis], win, poly, deriv=1, dt=dt)
        acc_obs[:, axis] = savgol_filter(xyz[:, axis], win, poly, deriv=2, dt=dt)

    # Compute gravitational acceleration (Earth + Moon)
    GM_earth = 398600.4418    # km³/s²
    GM_moon  = 4902.800066    # km³/s²

    acc_grav = np.zeros((n, 3))
    moon_xyz = np.zeros((n, 3))

    # JD for each record
    jd0 = 2451545.0  # J2000
    for i, r in enumerate(records):
        # JD from UTC
        delta = r['utc'] - datetime(2000, 1, 1, 12, 0, 0)
        jd = jd0 + delta.total_seconds() / 86400.0

        # Earth gravity
        r_vec = xyz[i]
        r_mag = np.linalg.norm(r_vec)
        acc_earth = -GM_earth * r_vec / r_mag**3

        # Moon gravity
        r_moon = moon_geocentric_xyz(jd)
        moon_xyz[i] = r_moon
        dr = xyz[i] - r_moon
        dr_mag = np.linalg.norm(dr)
        # Third-body acceleration = -GM * (dr/|dr|³ + r_moon/|r_moon|³)
        # (indirect term accounts for Moon's pull on Earth)
        acc_moon = -GM_moon * (dr / dr_mag**3 + r_moon / np.linalg.norm(r_moon)**3)

        acc_grav[i] = acc_earth + acc_moon

    # Residual = observed - gravitational
    acc_resid = acc_obs - acc_grav
    resid_mag = np.linalg.norm(acc_resid, axis=1)

    # Range rate (for reporting)
    ranges = np.array([r['range'] for r in records])
    range_rate = savgol_filter(ranges, win, poly, deriv=1, dt=dt)

    # Also compute a "raw" range acceleration for comparison
    range_accel_raw = np.diff(range_rate)  # forward differences

    # --- Detect burns: peaks in |residual acceleration| ---
    # Use a wide moving median as the "noise floor" and look for excursions
    hw = 30  # 5 hours each side
    margin = 8  # skip edge effects from S-G filter

    # Compute local median and MAD of residual magnitude
    from statistics import median
    resid_smooth = np.zeros(n)
    resid_mad = np.zeros(n)
    for i in range(n):
        lo = max(0, i - hw)
        hi = min(n, i + hw + 1)
        window = sorted(resid_mag[lo:hi])
        med = window[len(window) // 2]
        resid_smooth[i] = med
        devs = sorted(abs(v - med) for v in resid_mag[lo:hi])
        resid_mad[i] = devs[len(devs) // 2]

    # Normalised anomaly score
    anomaly = np.zeros(n)
    for i in range(n):
        if resid_mad[i] > 1e-12:
            anomaly[i] = (resid_mag[i] - resid_smooth[i]) / resid_mad[i]
        else:
            anomaly[i] = 0

    # Find peaks
    threshold = 6.0
    candidates = []
    for i in range(margin, n - margin):
        if anomaly[i] > threshold:
            # Local max within ±2
            local = anomaly[max(margin, i-2):min(n-margin, i+3)]
            if anomaly[i] >= max(local) * 0.90:
                candidates.append(i)

    # Cluster within 30 min
    events = []
    if candidates:
        current = [candidates[0]]
        for c in candidates[1:]:
            if c - current[-1] <= 3:
                current.append(c)
            else:
                events.append(current)
                current = [c]
        events.append(current)

    # --- Report ---
    print(f"\nRange: {records[0]['range']:.0f} → min {ranges.min():.0f} km → {records[-1]['range']:.0f} km")

    rmin_i = np.argmin(ranges)
    print(f"Closest approach: {records[rmin_i]['utc']} at {ranges[rmin_i]:.0f} km")

    # Moon closest approach
    sc_moon_dist = np.linalg.norm(xyz - moon_xyz, axis=1)
    moon_min_i = np.argmin(sc_moon_dist)
    print(f"Moon closest approach: {records[moon_min_i]['utc']} at {sc_moon_dist[moon_min_i]:.0f} km")

    print(f"\n{'='*80}")
    print(f"DETECTED {len(events)} CANDIDATE BURN EVENT(S)  (threshold = {threshold}σ above local median)")
    print(f"{'='*80}")

    for ei, event in enumerate(events):
        peak_idx = max(event, key=lambda i: anomaly[i])
        r = records[peak_idx]

        i_before = max(0, event[0] - 2)
        i_after = min(n - 1, event[-1] + 2)

        t0 = records[event[0]]['utc']
        t1 = records[event[-1]]['utc']
        dur = (t1 - t0).total_seconds() / 60

        # Residual acceleration components at peak
        ar = acc_resid[peak_idx]
        # Estimate ΔV by integrating residual acceleration over the event
        # (sum of residual × dt over the event window)
        dv = np.sum(acc_resid[event[0]:event[-1]+1], axis=0) * dt
        dv_mag = np.linalg.norm(dv)

        print(f"\n  BURN {ei+1}:")
        print(f"    Time:         {t0.strftime('%Y-%m-%d %H:%M')} UTC" +
              (f" – {t1.strftime('%H:%M')}" if dur > 0 else "") +
              f"  ({dur:.0f} min)")
        print(f"    Peak anomaly: {anomaly[peak_idx]:.1f}σ at {r['utc'].strftime('%H:%M')} UTC")
        print(f"    Range:        {r['range']:.0f} km   Moon dist: {sc_moon_dist[peak_idx]:.0f} km")
        print(f"    |Residual accel|: {resid_mag[peak_idx]*1000:.4f} m/s²")
        print(f"    Local background: {resid_smooth[peak_idx]*1000:.4f} m/s²")
        print(f"    Est. ΔV:      {dv_mag*1000:.1f} m/s  ({dv_mag:.4f} km/s)")
        print(f"    ΔV direction:  [{dv[0]:+.5f}, {dv[1]:+.5f}, {dv[2]:+.5f}] km/s (EME2000)")
        print(f"    Range rate:   {range_rate[i_before]:+.4f} → {range_rate[i_after]:+.4f} km/s")

    if not events:
        print("\n  No significant burns detected above threshold.")
        scored = [(i, anomaly[i]) for i in range(margin, n-margin)]
        scored.sort(key=lambda x: x[1], reverse=True)
        print(f"\n  Top 10 anomaly scores:")
        for idx, score in scored[:10]:
            print(f"    {records[idx]['utc']}  {score:.2f}σ  range={records[idx]['range']:.0f} km  |resid|={resid_mag[idx]*1000:.4f} m/s²")

    # Range extrema
    print(f"\n{'='*80}")
    print("RANGE EXTREMA")
    print(f"{'='*80}")
    for i in range(1, n - 1):
        if ranges[i] > ranges[i-1] and ranges[i] > ranges[i+1]:
            print(f"  APOGEE:  {records[i]['utc']}  {ranges[i]:>10.0f} km")
        elif ranges[i] < ranges[i-1] and ranges[i] < ranges[i+1]:
            print(f"  PERIGEE: {records[i]['utc']}  {ranges[i]:>10.0f} km")

    # --- Plot ---
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
    except ImportError:
        print("\nmatplotlib not available — skipping plots.")
        return

    times = [r['utc'] for r in records]

    fig, axes = plt.subplots(5, 1, figsize=(16, 16), sharex=True)
    fig.suptitle('Artemis II — Burn Detection\n(two-body gravity model subtracted, S-G differentiation)',
                 fontsize=13, fontweight='bold')

    # 1. Range + Moon distance
    ax = axes[0]
    ax.plot(times, ranges, 'b-', linewidth=0.8, label='Geocentric range')
    ax.plot(times, sc_moon_dist, 'gray', linewidth=0.8, alpha=0.6, label='Moon distance')
    ax.set_ylabel('km')
    ax.set_title('Range')
    ax.legend(fontsize=8)

    # 2. Range rate
    ax = axes[1]
    ax.plot(times, range_rate, 'b-', linewidth=0.8)
    ax.set_ylabel('km/s')
    ax.set_title('Range Rate')

    # 3. Residual acceleration magnitude
    ax = axes[2]
    ax.semilogy(times, resid_mag * 1000, 'r-', linewidth=0.6, alpha=0.7, label='|residual|')
    ax.semilogy(times, resid_smooth * 1000, 'k-', linewidth=1.2, alpha=0.5, label='local median')
    ax.set_ylabel('m/s²')
    ax.set_title('Non-gravitational Acceleration (observed − Earth − Moon gravity)')
    ax.legend(fontsize=8)

    # 4. Residual components
    ax = axes[3]
    ax.plot(times, acc_resid[:, 0]*1000, 'r-', linewidth=0.5, alpha=0.7, label='X')
    ax.plot(times, acc_resid[:, 1]*1000, 'g-', linewidth=0.5, alpha=0.7, label='Y')
    ax.plot(times, acc_resid[:, 2]*1000, 'b-', linewidth=0.5, alpha=0.7, label='Z')
    ax.set_ylabel('m/s² (EME2000)')
    ax.set_title('Residual Acceleration Components')
    ax.legend(fontsize=8)

    # 5. Anomaly score
    ax = axes[4]
    ax.plot(times, anomaly, 'k-', linewidth=0.8)
    ax.axhline(threshold, color='red', linestyle=':', alpha=0.7, label=f'threshold ({threshold}σ)')
    ax.set_ylabel('σ')
    ax.set_title('Anomaly Score (MAD-normalised)')
    ax.legend(fontsize=8)

    # Mark burns
    burn_colors = ['red', 'orange', 'darkviolet', 'brown', 'darkgreen', 'navy', 'deeppink']
    for ax in axes:
        ax.grid(True, alpha=0.3)
    for ei, event in enumerate(events):
        peak_idx = max(event, key=lambda i: anomaly[i])
        t = times[peak_idx]
        color = burn_colors[ei % len(burn_colors)]
        for ax in axes:
            ax.axvline(t, color=color, linestyle='--', alpha=0.8, linewidth=1.5)
        axes[0].annotate(f'Burn {ei+1}', xy=(t, axes[0].get_ylim()[1]*0.92),
                        fontsize=9, color=color, fontweight='bold', rotation=90,
                        va='top', ha='right')

    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter('%b %d\n%H:%M'))
    axes[-1].set_xlabel('UTC')
    plt.tight_layout()
    outpath = '/Users/karl/Dropbox/Software/Claude-Dev/ArtemisII/burn_detection.png'
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to {outpath}")
    plt.show()

if __name__ == '__main__':
    main()
