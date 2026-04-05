# Seestar / AAT / Horizons coordinate summary

## Core conclusions

- **Do not use geocentric coordinates** for a nearby spacecraft such as **Artemis II / Orion**.
- Use **topocentric** coordinates for the actual observing site and time.
- The remaining question is the **frame/epoch** expected by the telescope control system.

## Horizons output types

### Astrometric RA & DEC
- These are in an inertial reference frame such as **ICRF/J2000**.
- In Horizons, **astrometric** coordinates are **adjusted for light-time only**.
- That means they account for the fact that you see the spacecraft where it was when the light left it.
- They **do not** include the usual **stellar aberration due to the Earth's motion**.
- They also do **not** include the full **apparent-of-date** transformation.

### Apparent RA & DEC
- These are the coordinates appropriate to the **sky of date** as seen by the observer.
- They include:
  - precession and nutation to date
  - stellar aberration from observer motion
  - gravitational deflection
  - for Earth-based observing, the current-date equator/equinox frame

## Seestar conclusion

Empirically, Seestar appears to expect **J2000** custom coordinates and then convert them internally to **JNow**.

Therefore for Seestar custom entries, use:

**topocentric astrometric J2000 RA/Dec**

That avoids:
- using geocentric coordinates
- double application of precession/nutation
- feeding Seestar already-of-date coordinates

Do **not** enter Horizons **apparent RA/Dec** if Seestar is expecting **J2000** input.

## AAT conclusion

For a professional TCS such as the AAT, the correct input is likewise:

**topocentric astrometric J2000 RA/Dec**

Reason:
- the TCS expects **mean/catalog coordinates**
- it then performs its own internal conversion to the pointing frame of date
- that internal step is where precession/nutation and related pointing-frame updates are applied

Using Horizons **apparent** coordinates would risk doing the date-of-observation conversion twice.

Using Horizons **astrometric** coordinates does **not** double-count stellar aberration, because Horizons astrometric output does **not** include the observer-motion stellar aberration correction.

## Interpretation of Horizons wording

When Horizons says:

> "Astrometric RA & DEC ... Adjusted for light-time aberration only"

the correct interpretation is:

- this refers to the **finite light-travel time from the spacecraft**
- it accounts for the spacecraft's distance and motion during that light-travel time
- it does **not** mean the usual annual/stellar aberration due to the Earth's motion

So:

- **light-time correction** = yes
- **observer-motion stellar aberration** = no

## Practical recipe

For Artemis II / Orion in Horizons:

1. Use **Observer Table**
2. Use your **topocentric observing site**
3. Request **Astrometric RA & DEC**
4. Use the coordinates for the **exact minute** of observation
5. Enter those into Seestar or AAT if the system expects **J2000**

## One-line rule

If the TCS expects **J2000**, use **topocentric astrometric J2000** from Horizons, not **apparent-of-date** coordinates.
