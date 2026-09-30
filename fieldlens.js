/* ================================================================
   fieldlens.js -- exact paraxial electron optics of Glaser bell-shaped magnetic lenses.

   A round magnetic lens with on-axis field  B(z) = B0 / (1 + ((z - zc)/a)^2)  (Glaser's "Glockenfeld").
   In the frame that rotates with the electron (Larmor), the paraxial ray obeys

        r'' + g(z) r = 0,      g = k^2 / (a^2 (1 + zeta^2)^2),      zeta = (z - zc)/a,

   where k^2 = e B0^2 a^2 / (8 m phi*) is the dimensionless excitation (phi* the relativistic potential).
   The substitution  zeta = cot(phi),  r = u(phi)/sin(phi)  turns it into u'' + w^2 u = 0 with
   w = sqrt(1 + k^2), so the ray transfer matrix between ANY two planes is an exact closed form, and the
   Larmor rotation is  theta = k [atan(zeta2) - atan(zeta1)].

   Everything here is pure arithmetic on numbers: no DOM, no state. Lengths are in mm; a ray is (y, u = dy/dz).
   The same formulas were checked against numerical integration (scipy solve_ivp) to ~1e-11.
   ================================================================ */
"use strict";

const FieldLens = (function () {

  /** Solutions s1 = cos(w phi)/sin(phi), s2 = sin(w phi)/sin(phi) and their zeta-derivatives. */
  function basis(w, zeta) {
    const phi = Math.atan2(1, zeta), s = Math.sin(phi), c = Math.cos(phi);
    const cw = Math.cos(w * phi), sw = Math.sin(w * phi);
    return { s1: cw / s, s2: sw / s,
             d1: w * sw * s + cw * c,            // d s1 / d zeta
             d2: -w * cw * s + sw * c };         // d s2 / d zeta
  }

  /** Ray transfer matrix {A,B,C,D} from plane z1 to z2 through the field of a lens (a, zc, k2). */
  function matrix(k2, a, zc, z1, z2) {
    const w = Math.sqrt(1 + k2);
    const p = basis(w, (z1 - zc) / a), q = basis(w, (z2 - zc) / a);
    const f = -1 / w;                            // 1 / Wronskian, which is -w in zeta units
    const A = f * (q.s1 * p.d2 - q.s2 * p.d1);
    const Bz = f * (q.s2 * p.s1 - q.s1 * p.s2);
    const Cz = f * (q.d1 * p.d2 - q.d2 * p.d1);
    const Dz = f * (q.d2 * p.s1 - q.d1 * p.s2);
    return { A: A, B: a * Bz, C: Cz / a, D: Dz };    // zeta units -> mm
  }

  /** Larmor rotation of the meridional plane between z1 and z2, radians (sign is a convention). */
  function rotation(k2, a, zc, z1, z2) {
    return Math.sqrt(k2) * (Math.atan((z2 - zc) / a) - Math.atan((z1 - zc) / a));
  }

  /** Independent check: RK4 integration of r'' + g(z) r = 0, for any g(z) -- this is also the route for a
   *  tabulated field from a field solver. Returns {A,B,C,D} from z1 to z2. */
  function rk4Matrix(g, z1, z2, steps) {
    const h = (z2 - z1) / steps;
    const run = (y, u) => {
      let z = z1;
      const f = (zz, yy, uu) => [uu, -g(zz) * yy];
      for (let i = 0; i < steps; i++) {
        const k1 = f(z, y, u);
        const k2 = f(z + h / 2, y + h / 2 * k1[0], u + h / 2 * k1[1]);
        const k3 = f(z + h / 2, y + h / 2 * k2[0], u + h / 2 * k2[1]);
        const k4 = f(z + h, y + h * k3[0], u + h * k3[1]);
        y += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
        u += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
        z += h;
      }
      return [y, u];
    };
    const c1 = run(1, 0), c2 = run(0, 1);
    return { A: c1[0], B: c2[0], C: c1[1], D: c2[1] };
  }

  /** Decompose M = drift(dOut) . thin(f) . drift(dIn), exactly (always possible when C != 0). It is the
   *  cardinal-element form: a thin lens of focal length f located dIn past the entrance plane, with dOut
   *  left to the exit plane. dIn + dOut need NOT equal the length of the domain; the difference is how much
   *  a thick lens "shortens" the column, and it is what makes the thin-lens replacement exact. A lens that is
   *  switched off is a plain drift, for which any split will do: take the middle. */
  function cardinal(M, length) {
    if (Math.abs(M.C) < 1e-14) return { f: Infinity, dIn: length / 2, dOut: length / 2 };
    return { f: -1 / M.C, dIn: (M.D - 1) / M.C, dOut: (M.A - 1) / M.C };
  }

  /** Multiply 2x2 matrices, rightmost first:  compose(M2, M1) = M2 * M1. */
  function compose(M2, M1) {
    return { A: M2.A * M1.A + M2.B * M1.C, B: M2.A * M1.B + M2.B * M1.D,
             C: M2.C * M1.A + M2.D * M1.C, D: M2.C * M1.B + M2.D * M1.D };
  }
  const drift = d => ({ A: 1, B: d, C: 0, D: 1 });

  /* ---------------------------------------------------------------- lens objects */

  /** A lens owns the field between its entrance and exit planes [bm, bp]; beyond them the field is taken as
   *  zero. kind 'proj' is an ordinary round lens, whose strength is quoted as its projective focal length
   *  (object at infinity, whole domain). kind 'half' is the objective, whose specimen sits at the centre: its
   *  strength is quoted as the focal length of ONE HALF of the field (centre to exit), the published
   *  "pre-field f = post-field f". */
  function makeLens(spec) {
    const L = { id: spec.id, kind: spec.kind, zc: spec.zc, a: spec.a, bm: spec.bm, bp: spec.bp };
    L.matrix = (k2, z1, z2) => matrix(k2, L.a, L.zc, z1, z2);
    L.total = k2 => matrix(k2, L.a, L.zc, L.bm, L.bp);
    L.rotation = (k2, z1, z2) => rotation(k2, L.a, L.zc, z1, z2);
    // the strength the user sets, as a focal length in mm (Infinity when off)
    L.f = k2 => {
      const M = L.kind === 'half' ? matrix(k2, L.a, L.zc, L.bm, L.zc) : L.total(k2);
      return Math.abs(M.C) < 1e-16 ? Infinity : -1 / M.C;
    };
    // f falls from Infinity at k2 = 0 to a minimum and then rises again; the physical range is that first
    // branch. (Past it the solution is periodic and f dips again, so a blind search would find a later,
    // meaningless lobe.) Walk up from weak excitation until f stops falling, then refine that bracket.
    const fpos = k2 => { const v = L.f(k2); return v > 0 ? v : Infinity; };
    let lo = 0.02, hi = 0.04;
    while (fpos(hi) < fpos(lo) && hi < 40) { lo = hi; hi += 0.02; }
    lo = Math.max(1e-3, lo - 0.02);                    // bracket [lo, hi] now contains the minimum
    const g = 0.6180339887498949;
    let x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
    for (let i = 0; i < 100; i++) {
      if (fpos(x1) < fpos(x2)) { hi = x2; } else { lo = x1; }
      x1 = hi - g * (hi - lo); x2 = lo + g * (hi - lo);
    }
    L.k2Min = (lo + hi) / 2;
    L.fMin = L.f(L.k2Min);
    /** excitation k2 for a requested focal length; clamps to the strongest allowed (f >= fMin). */
    L.k2ForF = f => {
      if (!(f < 1e9)) return 0;
      if (f <= L.fMin) return L.k2Min;
      let a = 0, b = L.k2Min;
      for (let i = 0; i < 200; i++) { const m = (a + b) / 2; if (L.f(m) > f) a = m; else b = m; }
      return (a + b) / 2;
    };
    return L;
  }

  /** Choose the field width `a` so that the lens's strongest setting has exactly the focal length `fMin`
   *  (this also absorbs the small effect of cutting the field's 1/zeta^2 tails at the domain ends). */
  function calibrate(spec, fMinTarget) {
    let lo = 0.05, hi = 60;
    for (let i = 0; i < 80; i++) {
      const m = Math.sqrt(lo * hi);
      const fm = makeLens(Object.assign({}, spec, { a: m })).fMin;
      if (fm < fMinTarget) lo = m; else hi = m;
    }
    return Math.sqrt(lo * hi);
  }

  /** z values across [z1, z2], uniform in atan(zeta) so the dense part of the field gets the points. */
  function grid(L, z1, z2, n) {
    const t1 = Math.atan((z1 - L.zc) / L.a), t2 = Math.atan((z2 - L.zc) / L.a), out = new Array(n + 1);
    for (let i = 0; i <= n; i++) out[i] = L.zc + L.a * Math.tan(t1 + (t2 - t1) * i / n);
    out[0] = z1; out[n] = z2;
    return out;
  }

  /** The two basis trajectories through a lens, sampled on a grid, starting at zs[0]:
   *  a ray with initial (y0, u0) is  y0 * (y1, u1) + u0 * (y2, u2)  at each sample. */
  function basisCurves(L, k2, zs) {
    const n = zs.length, y1 = new Float64Array(n), u1 = new Float64Array(n),
          y2 = new Float64Array(n), u2 = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const M = i === 0 ? { A: 1, B: 0, C: 0, D: 1 } : matrix(k2, L.a, L.zc, zs[0], zs[i]);
      y1[i] = M.A; u1[i] = M.C; y2[i] = M.B; u2[i] = M.D;
    }
    return { z: zs, y1: y1, u1: u1, y2: y2, u2: u2 };
  }

  return { basis, matrix, rotation, rk4Matrix, cardinal, compose, drift, makeLens, calibrate, grid, basisCurves };
})();

if (typeof window !== 'undefined') window.FieldLens = FieldLens;
