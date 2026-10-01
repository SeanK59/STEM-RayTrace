/* ================================================================
   STEM column -- the single source of truth for both pages.

   index.html (the ray trace) loads this file. (The electron optics page keeps its own copy of the element
   positions in electron-column.js.)
   z in mm from the source crossover. f in mm. Ideal thin lenses throughout.
   ================================================================ */
"use strict";

// Element positions are MEASURED off the microscope. C3 -> PL1 is 610 mm and PL4 -> the
// EELS aperture is 400 mm with EL centred, both measured.
//
// The objective is the UltraSTEM 200's Riecke-Ruska condenser-objective: a 4 mm pole-piece gap
// and a 2 mm focal length for both pre- and post-field (Dellby et al., Eur. Phys. J. Appl. Phys.
// 54 (2011) 33505). It is modelled as two thin lenses on the pole faces, OL1 and OL2, 2 mm either
// side of the specimen -- so the specimen sits exactly on OL1's back focal plane and on OL2's
// front focal plane. A parallel beam into OL1 is focused on the specimen, and OL2 sends the
// specimen's image to infinity, with the diffraction pattern on its back focal plane at 789.
// (Krivanek et al., Ultramicroscopy 108 (2008) 179, gives 1.5 mm for the earlier 100 kV lens.)
// The VOA is 30 um in DIAMETER, the value Krivanek 2008 quotes as typical.
//
// C1 is off by default. C2 + C3 are the exact focused solve for 15 mrad from the bare point
// source: C3 collimates C2's crossover at 264, and alpha = beam radius at OL1 / f_OL1.
// PL1/PL2/PL4 are the exact solve for a 28.200 mm camera length at the EELS aperture with EL at
// a round 105 mm -- the gentlest four-lens family member, strongest lens 13.6 mm. At 39 mrad that
// camera length puts the direct beam at 1.10 mm radius on the aperture, which is the value
// measured on the real column. The projector digits are load-bearing: rounding them to 4 dp
// moves the EELS planes by ~1e-4 mm.
const COLUMN = [
  { id:'SRC', label:'Source',     z:    0, kind:'source' },
  { id:'C1',  label:'C1',         z:  120, kind:'lens',     f: Infinity },
  { id:'VOA', label:'VOA',        z:  180, kind:'aperture', r: 0.015 },
  { id:'C2',  label:'C2',         z:  240, kind:'lens',     f: 240 / 11 },
  { id:'C3',  label:'C3',         z:  300, kind:'lens',     f:  36 },
  { id:'OL1', label:'OL1',        z:  783, kind:'lens',     f:   2 },
  { id:'SMP', label:'Sample',     z:  785, kind:'plane',    marks:'image' },
  { id:'OL2', label:'OL2',        z:  787, kind:'lens',     f:   2 },
  { id:'PL1', label:'PL1',        z:  910, kind:'lens',     f:  56.0828398219 },
  { id:'PL2', label:'PL2',        z: 1030, kind:'lens',     f:  13.5835644681 },
  { id:'PL3', label:'PL3',        z: 1090, kind:'lens',     f: Infinity },
  { id:'PL4', label:'PL4',        z: 1150, kind:'lens',     f:  14.0309991241 },
  { id:'EL',  label:'EL',         z: 1350, kind:'lens',     f: 105 },
  { id:'EAP', label:'EELS ap.',   z: 1550, kind:'plane',    marks:'diff'  },
  { id:'EFP', label:'EELS focal', z: 1595, kind:'plane',    marks:'image' },
];

// Strongest lens allowed. Every round lens is f >= 4 mm; the objective, a far stronger
// immersion lens that normally runs at 2 mm, may go to 1.5 mm. These caps bound the solvers.
const P_MAX    = 1 / 4;       // round lenses: C1-C3, PL1-PL4, EL
const P_MAX_OL = 1 / 1.5;     // the objective's pre- and post-field
/** OL1 and OL2 are the two halves of ONE objective, always equally excited. */
const OBJECTIVE = ['OL1', 'OL2'];
const pMaxOf = id => OBJECTIVE.indexOf(id) >= 0 ? P_MAX_OL : P_MAX;

/** The nominal specimen height: midway between the pole faces. */
const COLUMN_Z_REF = (COLUMN.find(e => e.id === 'OL1').z + COLUMN.find(e => e.id === 'OL2').z) / 2;
/** Where the drawn column ends; the "parallel" test measures slope against this lever arm. */
const COLUMN_Z_END = 1655;

/** Paraxial transfer matrix [[A, B], [C, D]] from z1 to z2 (z2 >= z1), through every lens in
 *  (z1, z2] -- the same convention as the ray page's tracer, where a lens sitting exactly on the
 *  start plane has already acted and one on the end plane has just acted.
 *  power maps lens id -> 1/f in mm^-1 (0 = off). Returns {A, B, C, D}. */
function transfer(power, z1, z2) {
  let A = 1, B = 0, C = 0, D = 1, z = z1;
  for (const el of COLUMN) {
    if (el.kind !== 'lens' || !(el.z > z1 + 1e-9 && el.z <= z2 + 1e-9)) continue;
    const d = el.z - z;                    // drift, then the lens
    A += d * C; B += d * D; z = el.z;
    const p = power[el.id] || 0;
    C -= p * A; D -= p * B;
  }
  const d = z2 - z;
  return { A: A + d * C, B: B + d * D, C: C, D: D };
}

/** What the condenser delivers to the specimen, from the lens powers alone.
 *
 *  The illumination is a point source clipped by the VOA, so every ray is a multiple of the
 *  unit-slope axial ray from z = 0, and u0max -- the slope of the ray that grazes the VOA edge --
 *  scales it. alpha is that marginal ray's slope inside the objective gap (no lens sits between
 *  OL1 and OL2, so the slope is the same everywhere in it); under parallel illumination it is the
 *  angle OL2 focuses the beam to instead, with the illuminated radius alpha * f_OL2. The ray page's
 *  analyse() defines them identically, and its self-test asserts the two agree.
 *  Returns { u0max, alpha (rad), collimated, spotR (mm, beam radius at the reference plane),
 *            zFocus (mm, where the axial ray crosses the axis if it did so in the gap's straight
 *            line; Infinity when parallel) }, or null if the VOA sits on a crossover. */
function illumination(power, voaR) {
  const zVOA = COLUMN.find(e => e.id === 'VOA').z;
  const yVOA = transfer(power, 0, zVOA).B;             // the axial ray (0, 1) arrives at height B
  if (!(Math.abs(yVOA) > 1e-9)) return null;
  const u0max = voaR / Math.abs(yVOA);
  const ref = transfer(power, 0, COLUMN_Z_REF);
  const y = ref.B, u = ref.D;
  const collimated = Math.abs(u) * (COLUMN_Z_END - COLUMN_Z_REF) < 1e-3 * Math.abs(y);
  const pOL2 = power.OL2 || 0;
  const alpha = (collimated && pOL2 > 0) ? Math.abs(y) * u0max * pOL2 : Math.abs(u) * u0max;
  return { u0max: u0max, alpha: alpha, collimated: collimated, spotR: Math.abs(y) * u0max,
           zFocus: u !== 0 ? COLUMN_Z_REF - y / u : Infinity };
}
