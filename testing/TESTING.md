# Artemis II Tracker -- Testing & Validation

This folder contains scripts, test data, and output plots used to validate
the topocentric coordinate calculations in `artemis2_track.html`.

## Summary of findings

1. **JS topocentric code vs Astropy**: The app's precession/nutation/GAST
   transformation chain (Lieske 1976 + IAU 1980 4-term nutation) agrees with
   Astropy's full GCRS/ITRS machinery to **~0.2 arcsec RMS** across the full
   mission. Verified in both northern and southern hemispheres -- no
   hemisphere-dependent bugs.

2. **OEM vs Horizons trajectory divergence**: Comparing OEM-derived positions
   against live JPL Horizons shows ~2-6 arcsec agreement in mid-mission
   (Apr 3-8), with exponential divergence near re-entry (Apr 10-11). This is
   expected -- different trajectory nav solutions (OEM generated days earlier)
   predict increasingly different paths as re-entry approaches.

## Scripts

### `validate_topo_vs_astropy.py`

Validates the JS topocentric calculation against Astropy at a single observer
(Melbourne). Reads geocentric EME2000 XYZ from the OEM data embedded in
`artemis2_track.html`, computes topocentric RA/Dec two ways:

1. Python replica of the JS code (Lieske precession + IAU 1980 nutation + GAST)
2. Astropy GCRS -> ITRS -> observer subtraction

Prints a table of JS vs Astropy RA/Dec at sampled times, plus RMS/max
separation statistics.

**Requirements:** `numpy`, `astropy`

```bash
cd testing
python3 validate_topo_vs_astropy.py
```

**Expected output:** Max separation < 1 arcsec, RMS ~ 0.2 arcsec.

---

### `validate_topo_hemisphere.py`

Same validation as above, but run at **two** observer locations:
- Siding Spring Observatory, Australia (31.3°S, 149.1°E) -- southern hemisphere
- UK (53°N, 1°E) -- northern hemisphere

Confirms there is no hemisphere-dependent bug in the WGS84 geodetic-to-ECEF
conversion or the GAST Earth rotation.

**Requirements:** `numpy`, `astropy`

```bash
cd testing
python3 validate_topo_hemisphere.py
```

**Expected output:** Both hemispheres agree with Astropy to < 1 arcsec RMS.

---

### `compare_oem_vs_horizons.py`

Compares the app's OEM-derived topocentric output against JPL Horizons
topocentric output at three observer locations. The OEM data is interpolated
onto the Horizons 10-minute time grid before differencing, to avoid artefacts
from timestamp misalignment.

Generates three 3-panel plots (dRA, dDec, total separation vs time) in
`plots/`.

**Requirements:** `numpy`, `matplotlib`

```bash
cd testing
python3 compare_oem_vs_horizons.py
```

**Expected output:** Mid-mission agreement of ~2-6 arcsec (trajectory vintage
differences), with exponential divergence near re-entry. Plots saved to
`plots/diff_uk.png`, `plots/diff_melbourne.png`, `plots/diff_sso.png`.

---

## Test data files (`data/`)

All files are "Save as Text" exports from `artemis2_track.html` in topocentric
ICRF RA/Dec format (degrees). Each file records the data source, generation
timestamp, and observer location in its header.

### OEM-derived (computed by the app from embedded CCSDS OEM ephemeris)

| File | OEM source | Generated | Observer |
|------|-----------|-----------|----------|
| `oem_v3_melbourne.txt` | Artemis_II_OEM_2026_04_02_to_EI_v3.asc | 2026-04-02 | Melbourne (37.9°S, 145.0°E) |
| `oem_v3_uk.txt` | Artemis_II_OEM_2026_04_02_to_EI_v3.asc | 2026-04-02 | UK (53°N, 1°E) |
| `oem_apr04_melbourne.txt` | NASA_Support_Orion_0404_OSA_OEM_NO_COV.asc | 2026-04-04 | Melbourne (37.9°S, 145.0°E) |
| `oem_apr04_sso.txt` | NASA_Support_Orion_0404_OSA_OEM_NO_COV.asc | 2026-04-04 | Siding Spring (31.3°S, 149.1°E) |
| `oem_apr04_uk.txt` | NASA_Support_Orion_0404_OSA_OEM_NO_COV.asc | 2026-04-04 | UK (53°N, 1°E) |

### Horizons-derived (fetched live from JPL Horizons via the app)

| File | Generated | Observer |
|------|-----------|----------|
| `horizons_melbourne.txt` | 2026-04-05 | Melbourne (37.9°S, 145.0°E) |
| `horizons_sso.txt` | 2026-04-05 | Siding Spring (31.3°S, 149.1°E) |
| `horizons_uk.txt` | 2026-04-05 | UK (53°N, 1°E) |

## Output plots (`plots/`)

| Plot | Description |
|------|-------------|
| `diff_uk.png` | OEM v3 (Apr 2) vs Horizons (Apr 5), UK observer |
| `diff_melbourne.png` | OEM v3 (Apr 2) vs Horizons (Apr 5), Melbourne observer |
| `diff_sso.png` | OEM Apr-04 vs Horizons (Apr 5), Siding Spring observer |

Each plot has three panels:
1. **dRA*cos(dec)** in arcseconds (blue) -- RA difference projected on sky
2. **dDec** in arcseconds (red) -- declination difference
3. **Total separation** in arcseconds on log scale (black)
