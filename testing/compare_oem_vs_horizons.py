#!/usr/bin/env python3
"""
Compare OEM-computed topocentric RA/Dec vs JPL Horizons topocentric RA/Dec.

Reads two app-export text files (OEM-derived and Horizons-derived),
interpolates the OEM data onto the Horizons time grid, computes angular
differences, and generates a 3-panel difference-vs-time plot for each
observer location.

Requires: numpy, matplotlib

Run from the testing/ directory:
    python3 compare_oem_vs_horizons.py
"""

import os
import re
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta

def parse_app_output(path):
    """Parse app output in DEG format. Returns list of (datetime, ra, dec), plus header info."""
    records = []
    source = observer = generated = ''
    with open(path) as f:
        for line in f:
            if line.startswith('Source:'):
                source = line.strip()[7:].strip()
            if 'Generated:' in line:
                m = re.search(r'Generated:\s*(\S+)', line)
                if m: generated = m.group(1)
            if 'Observer:' in line:
                m = re.search(r'Observer:\s*lat=([\d.+-]+)\s*lon=([\d.+-]+)', line)
                if m: observer = f"{m.group(1)}°, {m.group(2)}°"
            parts = line.strip().split()
            utc = None
            for p in parts:
                if p.startswith('2026-'):
                    utc = p; break
            if not utc: continue
            idx = parts.index(utc)
            try:
                ra = float(parts[idx+1])
                dec = float(parts[idx+2])
            except (IndexError, ValueError):
                continue
            fmt = '%Y-%m-%dT%H:%M:%S' if len(utc) >= 19 else '%Y-%m-%dT%H:%M'
            dt = datetime.strptime(utc[:min(len(utc), 19)], fmt)
            if '.' in utc and len(utc) > 19:
                frac = float('0' + utc[19:])
                dt = dt + timedelta(seconds=frac)
            records.append((dt, ra, dec))
    records.sort(key=lambda r: r[0])
    return records, source, observer, generated


def interpolate_radec(records, query_times):
    """Interpolate RA/Dec at query_times. Handles RA wraparound."""
    ref_epoch = records[0][0]
    t_sec = np.array([(r[0] - ref_epoch).total_seconds() for r in records])
    ra_arr = np.array([r[1] for r in records])
    dec_arr = np.array([r[2] for r in records])
    ra_unwrap = np.unwrap(np.radians(ra_arr))

    results = []
    for qt in query_times:
        qs = (qt - ref_epoch).total_seconds()
        if qs < t_sec[0] or qs > t_sec[-1]:
            results.append(None)
            continue
        ra_interp = np.interp(qs, t_sec, ra_unwrap)
        dec_interp = np.interp(qs, t_sec, dec_arr)
        results.append((math.degrees(ra_interp) % 360, dec_interp))
    return results


def compare_and_plot(file_oem, file_hz, location_name, outpath):
    """Run interpolated comparison and generate a 3-panel difference plot."""
    recs_oem, src_oem, obs_oem, gen_oem = parse_app_output(file_oem)
    recs_hz, src_hz, obs_hz, gen_hz = parse_app_output(file_hz)

    print(f"\n{'='*70}")
    print(f"  {location_name}")
    print(f"  OEM:      {src_oem}  (generated {gen_oem})")
    print(f"  Horizons: {src_hz}  (generated {gen_hz})")
    print(f"  Records:  OEM={len(recs_oem)}  Horizons={len(recs_hz)}")
    print(f"{'='*70}")

    # Interpolate OEM onto the Horizons time grid
    grid_times = [r[0] for r in recs_hz]
    grid_radec = [(r[1], r[2]) for r in recs_hz]
    interp_oem = interpolate_radec(recs_oem, grid_times)

    times, dra_list, ddec_list, dtot_list = [], [], [], []
    for i, qt in enumerate(grid_times):
        if interp_oem[i] is None:
            continue
        a_ra, a_dec = interp_oem[i]
        b_ra, b_dec = grid_radec[i]

        dra_raw = a_ra - b_ra
        if dra_raw > 180: dra_raw -= 360
        if dra_raw < -180: dra_raw += 360
        mean_dec = (a_dec + b_dec) / 2
        dra_sky = dra_raw * 3600 * math.cos(math.radians(mean_dec))
        ddec = (a_dec - b_dec) * 3600
        dtot = math.sqrt(dra_sky**2 + ddec**2)

        times.append(qt)
        dra_list.append(dra_sky)
        ddec_list.append(ddec)
        dtot_list.append(dtot)

    rms_dra = math.sqrt(sum(x**2 for x in dra_list) / len(dra_list))
    rms_ddec = math.sqrt(sum(x**2 for x in ddec_list) / len(ddec_list))
    rms_tot = math.sqrt(sum(x**2 for x in dtot_list) / len(dtot_list))

    print(f"  Compared: {len(times)} points")
    print(f"  dRA*cos(dec):  min={min(dra_list):+.1f}\"  max={max(dra_list):+.1f}\"  rms={rms_dra:.1f}\"")
    print(f"  dDec:          min={min(ddec_list):+.1f}\"  max={max(ddec_list):+.1f}\"  rms={rms_ddec:.1f}\"")
    print(f"  Total sep:     min={min(dtot_list):.1f}\"   max={max(dtot_list):.1f}\"   rms={rms_tot:.1f}\"")

    # Worst 5
    ranked = sorted(zip(dtot_list, times, dra_list, ddec_list), reverse=True)
    print(f"  Worst 5:")
    for tot, t, dra, ddec in ranked[:5]:
        print(f"    {t.strftime('%Y-%m-%d %H:%M')}  dRA={dra:+.1f}\"  dDec={ddec:+.1f}\"  total={tot:.1f}\"")

    # ── Plot ──
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)

    axes[0].plot(times, dra_list, 'b-', linewidth=0.7)
    axes[0].set_ylabel('dRA*cos(dec) (arcsec)')
    axes[0].axhline(0, color='gray', linewidth=0.5)
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(times, ddec_list, 'r-', linewidth=0.7)
    axes[1].set_ylabel('dDec (arcsec)')
    axes[1].axhline(0, color='gray', linewidth=0.5)
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(times, dtot_list, 'k-', linewidth=0.7)
    axes[2].set_ylabel('Total separation (arcsec)')
    axes[2].set_yscale('log')
    axes[2].grid(True, alpha=0.3)
    axes[2].set_xlabel('UTC')

    for ax in axes:
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d\n%H:%M'))
        ax.xaxis.set_major_locator(mdates.DayLocator())

    oem_short = src_oem.replace('.asc', '').replace('Artemis_II_OEM_2026_04_02_to_EI_', 'OEM ').replace('NASA_Support_Orion_0404_OSA_OEM_NO_COV', 'OEM Apr-04')
    fig.suptitle(f'OEM vs JPL Horizons — {location_name}\n'
                 f'{oem_short} ({gen_oem[:10]}) vs Horizons ({gen_hz[:10]})',
                 fontsize=13)
    fig.autofmt_xdate()
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  Plot → {outpath}")


# ═══════════════════════════════════════════════════════════════
# Run all three comparisons
# ═══════════════════════════════════════════════════════════════

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
data = os.path.join(SCRIPT_DIR, 'data')
plots = os.path.join(SCRIPT_DIR, 'plots')
os.makedirs(plots, exist_ok=True)

compare_and_plot(
    os.path.join(data, 'oem_v3_uk.txt'),
    os.path.join(data, 'horizons_uk.txt'),
    'UK (53.0°N, 1.0°E)',
    os.path.join(plots, 'diff_uk.png'))

compare_and_plot(
    os.path.join(data, 'oem_v3_melbourne.txt'),
    os.path.join(data, 'horizons_melbourne.txt'),
    'Melbourne (37.9°S, 145.0°E)',
    os.path.join(plots, 'diff_melbourne.png'))

compare_and_plot(
    os.path.join(data, 'oem_apr04_sso.txt'),
    os.path.join(data, 'horizons_sso.txt'),
    'Siding Spring (31.3°S, 149.1°E)',
    os.path.join(plots, 'diff_sso.png'))
