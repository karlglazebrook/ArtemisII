# Artemis II Mission Tracker

A self-contained, single-file browser application for tracking NASA's Artemis II
crewed lunar mission in real time. Displays the spacecraft's position on the sky
as seen from any observer location on Earth.

![Artemis II Tracker](https://img.shields.io/badge/mission-Artemis_II-blue)
![Self Contained](https://img.shields.io/badge/dependencies-none-green)

## Features

- **Single HTML file** -- no server, no build step, no dependencies. Just open
  `artemis2_track.html` in a browser.
- **Topocentric ICRF RA/Dec** -- observer-dependent coordinates computed from
  geocentric EME2000 (J2000) state vectors using Lieske precession, IAU 1980
  nutation, and GAST Earth rotation.
- **Dual data sources**:
  - **Embedded CCSDS OEM** -- NASA trajectory ephemeris baked into the file for
    offline use.
  - **JPL Horizons live fetch** -- queries Horizons API for the latest
    trajectory solution (NAIF ID -1024), with CORS proxy fallback.
- **Configurable observer** -- set latitude, longitude, altitude, or use
  browser geolocation. Defaults to Siding Spring Observatory.
- **Multiple views**:
  - Data table with RA/Dec in degrees and HMS/DMS formats
  - Sky track plot (RA/Dec over time)
  - Real-time sky view (alt/az polar plot with horizon)
  - Panoramic horizon view
  - Range and angular rate charts
- **Text export** with full provenance (data source, generation time, observer).

## Quick start

```bash
open artemis2_track.html
# or
python3 -m http.server 8765   # then visit http://localhost:8765/artemis2_track.html
```

For live Horizons data during development, disable the browser's same-origin
policy:

```bash
# Chrome
open -na "Google Chrome" --args --disable-web-security --user-data-dir=/tmp/chrome-cors-dev

# Safari
# Develop menu → Disable Cross-Origin Restrictions
```

## Files

| File | Description |
|------|-------------|
| `artemis2_track.html` | The complete tracker application (single file, ~650KB) |
| `Artemis_II_OEM_2026_04_02_to_EI_v3.asc` | CCSDS OEM ephemeris, nav solution from Apr 2 |
| `parse_oem.js` | Standalone OEM parser (Node.js, used during development) |
| `detect_burns.py` | Script to detect propulsive manoeuvres from OEM velocity data |
| `proxy_server.py` | Local CORS proxy for Horizons API (development aid) |
| `testing/` | Validation scripts, test data, and comparison plots |

## Topocentric calculation

The app converts geocentric EME2000 state vectors to observer-centric RA/Dec
via:

1. **WGS84 geodetic to ECEF** -- observer lat/lon/alt to Earth-fixed XYZ
2. **GAST rotation** -- ECEF to true-of-date equatorial (Earth rotation angle
   including equation of the equinoxes)
3. **Inverse nutation** (IAU 1980, 4-term) -- true-of-date to mean-of-date
4. **Inverse precession** (Lieske 1976) -- mean-of-date to J2000/EME2000
5. **Topocentric vector** -- subtract observer position from spacecraft
   position, convert to RA/Dec

This pipeline has been validated against Astropy's full GCRS/ITRS
transformation to **~0.2 arcsec RMS** across the full mission duration and in
both hemispheres. See [`testing/TESTING.md`](testing/TESTING.md) for details.

## Testing

The `testing/` folder contains:

- **`validate_topo_vs_astropy.py`** -- validates the JS topocentric code
  against Astropy (requires `numpy`, `astropy`)
- **`validate_topo_hemisphere.py`** -- confirms no hemisphere-dependent bugs
- **`compare_oem_vs_horizons.py`** -- compares OEM vs Horizons trajectories
  at three observer locations with time interpolation (requires `numpy`,
  `matplotlib`)
- **`data/`** -- app text exports used as test inputs
- **`plots/`** -- generated difference plots

```bash
cd testing
python3 validate_topo_vs_astropy.py
python3 validate_topo_hemisphere.py
python3 compare_oem_vs_horizons.py
```

See [`testing/TESTING.md`](testing/TESTING.md) for full documentation.

## Licence

This project is provided as-is for educational and amateur astronomy purposes.
Trajectory data is sourced from NASA/JPL.
