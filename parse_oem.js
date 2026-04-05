#!/usr/bin/env node
/**
 * parse_oem.js
 *
 * Reads a CCSDS OEM ephemeris file (REF_FRAME = EME2000) and outputs
 * geocentric RA / DEC in J2000 (epoch & equinox 2000.0) vs UTC.
 *
 * Because the OEM states are already in the EME2000 frame
 * (Earth Mean Equator & Equinox of J2000), the conversion is purely
 * a rectangular-to-spherical coordinate transformation — no library needed.
 *
 *   RA  = atan2(Y, X)            [0, 360) deg  or [0h, 24h)
 *   Dec = asin(Z / |r|)          (−90, +90) deg
 *
 * Usage:
 *   node parse_oem.js <ephemeris.asc> [--csv]
 *
 * Output (default JSON):
 *   [{ utc, jd, x_km, y_km, z_km, range_km, ra_deg, dec_deg, ra_hms, dec_dms }, ...]
 *
 * Output (--csv):
 *   UTC, JD, X_km, Y_km, Z_km, Range_km, RA_deg, Dec_deg, RA_hms, Dec_dms
 */

'use strict';

const fs   = require('fs');
const path = require('path');

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Convert degrees to radians */
const DEG = Math.PI / 180;

/** Julian Date of J2000.0 (2000-Jan-01 12:00 TT ≈ UTC) */
const JD_J2000 = 2451545.0;

/**
 * Parse an ISO 8601 UTC timestamp from the OEM
 * e.g. "2026-04-02T03:07:49.583"
 * Returns a JavaScript Date object (UTC).
 */
function parseOEMDate(s) {
  return new Date(s + 'Z');
}

/**
 * Compute Julian Date from a JavaScript Date (UTC).
 */
function dateToJD(d) {
  return d.getTime() / 86400000.0 + 2440587.5;
}

/**
 * Format decimal degrees as ±DD°MM'SS.ss"
 */
function degToDMS(deg) {
  const sign = deg < 0 ? '-' : '+';
  const abs  = Math.abs(deg);
  const d    = Math.floor(abs);
  const mf   = (abs - d) * 60;
  const m    = Math.floor(mf);
  const s    = (mf - m) * 60;
  return `${sign}${String(d).padStart(2,'0')}°${String(m).padStart(2,'0')}'${s.toFixed(2).padStart(5,'0')}"`;
}

/**
 * Format decimal degrees of RA as HH:MM:SS.ss
 */
function degToHMS(deg) {
  const h_total = ((deg % 360) + 360) % 360 / 15;
  const h  = Math.floor(h_total);
  const mf = (h_total - h) * 60;
  const m  = Math.floor(mf);
  const s  = (mf - m) * 60;
  return `${String(h).padStart(2,'0')}h${String(m).padStart(2,'0')}m${s.toFixed(2).padStart(5,'0')}s`;
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

function parseOEM(text) {
  const lines   = text.split('\n');
  const records = [];
  let inMeta    = false;
  let inData    = false;
  let frame     = null;

  for (const raw of lines) {
    const line = raw.trim();

    if (line.startsWith('META_START')) { inMeta = true;  inData = false; continue; }
    if (line.startsWith('META_STOP'))  { inMeta = false; inData = true;  continue; }

    if (inMeta) {
      const m = line.match(/^REF_FRAME\s*=\s*(\S+)/);
      if (m) frame = m[1];
      continue;
    }

    if (!inData) continue;
    if (!line || line.startsWith('COMMENT') || line.startsWith('#')) continue;

    // Data line: YYYY-MM-DDTHH:MM:SS.sss  X Y Z Vx Vy Vz
    const parts = line.split(/\s+/);
    if (parts.length < 7) continue;
    const dateStr = parts[0];
    if (!dateStr.match(/^\d{4}-\d{2}-\d{2}T/)) continue;

    records.push({
      dateStr,
      x: parseFloat(parts[1]),
      y: parseFloat(parts[2]),
      z: parseFloat(parts[3]),
      vx: parseFloat(parts[4]),
      vy: parseFloat(parts[5]),
      vz: parseFloat(parts[6]),
    });
  }

  if (frame && frame !== 'EME2000') {
    process.stderr.write(`WARNING: REF_FRAME is "${frame}", not EME2000. ` +
      `RA/Dec values will be incorrect without a frame transformation.\n`);
  }

  return records;
}

// ---------------------------------------------------------------------------
// Coordinate conversion: EME2000 XYZ → geocentric RA, Dec (J2000)
// ---------------------------------------------------------------------------

function xyzToRADec(x, y, z) {
  const r   = Math.sqrt(x*x + y*y + z*z);          // km
  const ra  = Math.atan2(y, x) / DEG;               // deg, −180..+180
  const dec = Math.asin(z / r) / DEG;               // deg

  return {
    range_km : r,
    ra_deg   : ((ra % 360) + 360) % 360,            // 0..360
    dec_deg  : dec,
  };
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

function main() {
  const args    = process.argv.slice(2);
  const csvMode = args.includes('--csv');
  const file    = args.find(a => !a.startsWith('--'));

  if (!file) {
    process.stderr.write(`Usage: node parse_oem.js <file.asc> [--csv]\n`);
    process.exit(1);
  }

  const text    = fs.readFileSync(path.resolve(file), 'utf8');
  const records = parseOEM(text);

  const results = records.map(r => {
    const date  = parseOEMDate(r.dateStr);
    const jd    = dateToJD(date);
    const coord = xyzToRADec(r.x, r.y, r.z);

    return {
      utc      : r.dateStr,
      jd       : jd,
      x_km     : r.x,
      y_km     : r.y,
      z_km     : r.z,
      range_km : coord.range_km,
      ra_deg   : coord.ra_deg,
      dec_deg  : coord.dec_deg,
      ra_hms   : degToHMS(coord.ra_deg),
      dec_dms  : degToDMS(coord.dec_deg),
    };
  });

  if (csvMode) {
    const header = 'UTC,JD,X_km,Y_km,Z_km,Range_km,RA_deg,Dec_deg,RA_hms,Dec_dms';
    process.stdout.write(header + '\n');
    for (const r of results) {
      process.stdout.write(
        `${r.utc},${r.jd.toFixed(6)},` +
        `${r.x_km.toFixed(3)},${r.y_km.toFixed(3)},${r.z_km.toFixed(3)},` +
        `${r.range_km.toFixed(3)},${r.ra_deg.toFixed(6)},${r.dec_deg.toFixed(6)},` +
        `"${r.ra_hms}","${r.dec_dms}"\n`
      );
    }
  } else {
    process.stdout.write(JSON.stringify(results, null, 2) + '\n');
  }

  process.stderr.write(`Processed ${results.length} records\n`);
  process.stderr.write(`Time span: ${results[0].utc} → ${results[results.length-1].utc}\n`);
  process.stderr.write(`RA range: ${Math.min(...results.map(r=>r.ra_deg)).toFixed(2)}° – ` +
    `${Math.max(...results.map(r=>r.ra_deg)).toFixed(2)}°\n`);
  process.stderr.write(`Dec range: ${Math.min(...results.map(r=>r.dec_deg)).toFixed(2)}° – ` +
    `${Math.max(...results.map(r=>r.dec_deg)).toFixed(2)}°\n`);
}

main();
