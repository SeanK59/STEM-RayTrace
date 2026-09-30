/* ================================================================
   electron-column.js -- the STEM column as real lens FIELDS (Glaser bell fields), not thin lenses.

   Depends on fieldlens.js. Pure logic, no DOM. z in mm from the source; a ray is (y, u = dy/dz) in the frame
   that rotates with the electron.

   Element positions are the same measured values the thin-lens page uses, so the two pages describe one
   column. Each round lens owns the field between its entrance and exit planes (cut at the midpoint to its
   neighbours, or at a fixed cap where they are far away); beyond that the field is taken as zero. The objective
   is ONE lens, centred on the specimen plane, whose strength is quoted as the focal length of one half of its
   field (the published pre-field = post-field f).

   What the page needs from here:
     - the exact thick-lens optics: transfer matrices between any two planes, sampled trajectories,
       Larmor rotation, plane finding                                               (exact*, below)
     - an EFFECTIVE thin-lens column: every thick lens is, exactly, a thin lens of focal length f at a
       principal-plane position, with the drifts shortened or lengthened to match. The page's validated
       thin-lens solvers run on that column, and the answer is polished against the exact matrices.
   ================================================================ */
"use strict";

const ElectronColumn = (function () {
  const FL = FieldLens;

  /* ------------------------------------------------------------ the physical column */
  const Z_MIN = -60, Z_MAX = 1655;
  const Z_SRC = 0, Z_VOA = 180, Z_OL = 785, Z_EAP = 1550, Z_EFP = 1595;
  const VOA_R = 0.015;                                   // mm: the VOA is 30 um in diameter

  // Lens field centres, in order. 'proj' lenses are quoted by their projective focal length; 'half' (the
  // objective) by the focal length of one half of its field.
  const LENSES = [
    { id: 'C1',  zc:  120, kind: 'proj' },
    { id: 'C2',  zc:  240, kind: 'proj' },
    { id: 'C3',  zc:  300, kind: 'proj' },
    { id: 'OL',  zc:  785, kind: 'half' },
    { id: 'PL1', zc:  910, kind: 'proj' },
    { id: 'PL2', zc: 1030, kind: 'proj' },
    { id: 'PL3', zc: 1090, kind: 'proj' },
    { id: 'PL4', zc: 1150, kind: 'proj' },
    { id: 'EL',  zc: 1350, kind: 'proj' },
  ];
  // the strongest setting of each lens. The objective's 2 mm is the UltraSTEM 200's published f (Dellby et
  // al. 2011); 4 mm is the smallest focal length of the round lenses' pole pieces (Krivanek et al. 2008).
  const F_MIN = { proj: 4, half: 2 };
  const CAP = { proj: 40, half: 30 };                    // mm: fields are cut this far from a lens at most

  const order = LENSES.map(l => l.id);

  // field domains: halfway to the neighbouring lens, or the cap
  LENSES.forEach((l, i) => {
    const prev = i > 0 ? LENSES[i - 1].zc : -Infinity, next = i < LENSES.length - 1 ? LENSES[i + 1].zc : Infinity;
    l.bm = Math.max(l.zc - CAP[l.kind], (l.zc + prev) / 2);
    l.bp = Math.min(l.zc + CAP[l.kind], (l.zc + next) / 2);
  });

  // Field widths. `a` is chosen so each lens's strongest setting has exactly F_MIN: this pins the
  // published focal lengths and absorbs the tiny effect of cutting the 1/zeta^2 tails at the domain ends. The
  // objective comes out at a = 2.013 mm, i.e. FWHM 4.03 mm against its 4 mm pole gap. The selftest recomputes
  // these from scratch and checks them.
  const WIDTH = {};
  const build = () => {
    const out = {};
    for (const l of LENSES) {
      const spec = { id: l.id, kind: l.kind, zc: l.zc, bm: l.bm, bp: l.bp };
      const a = WIDTH[l.id] !== undefined ? WIDTH[l.id] : FL.calibrate(spec, F_MIN[l.kind]);
      WIDTH[l.id] = a;
      out[l.id] = FL.makeLens(Object.assign({}, spec, { a: a }));
    }
    return out;
  };
  const lens = build();

  /* ------------------------------------------------------------ the specimen's nominal plane */
  // The strongest objective focuses a PARALLEL beam not on the centre of its field but a few um above it (its
  // pre-field principal plane sits 2.0057 mm ahead of the centre for f = 2.0000). That is where the real
  // microscope runs its probe -- Dellby et al. 2011: "the sample-level beam crossover a few um above the center
  // of the objective lens" -- and it is where the specimen's nominal plane belongs: measured from there, a
  // parallel beam focuses exactly on it, instead of needing a slightly diverging beam (which the 30 um-wide beam
  // makes hypersensitive). Found by bisection: the height of the unit-slope... of the parallel ray vanishes there.
  const Z_REF = (function () {
    const L = lens.OL, k2 = L.k2Min;
    let lo = L.zc - 0.5, hi = L.zc + 0.5;
    const A = z => L.matrix(k2, L.bm, z).A;           // height at z of a ray entering parallel with unit height
    for (let i = 0; i < 100; i++) { const m = (lo + hi) / 2; if (A(lo) * A(m) <= 0) hi = m; else lo = m; }
    return (lo + hi) / 2;
  })();

  /* ------------------------------------------------------------ excitation <-> strength */
  // f -> k2 is a bisection, so cache it: a slider drag revisits the same few values constantly.
  const k2cache = {};
  function k2Of(id, power) {
    if (!(power > 0)) return 0;
    const c = k2cache[id];
    if (c && c.p === power) return c.k2;
    const k2 = lens[id].k2ForF(1 / power);
    k2cache[id] = { p: power, k2: k2 };
    return k2;
  }
  /** strongest allowed power of a lens (1 / f_min) */
  const pMax = id => 1 / lens[id].fMin;

  /* ------------------------------------------------------------ exact optics */
  const I = { A: 1, B: 0, C: 0, D: 1 };
  const mul = FL.compose;

  /** Exact transfer matrix from z1 to z2 (z2 >= z1) through every field in between, for `powers`
   *  ({id: 1/f}). A lens that is off, or a stretch with no field, is a plain drift. */
  function chain(powers, z1, z2) {
    let M = I, z = z1;
    for (const id of order) {
      const L = lens[id], lo = Math.max(z, L.bm), hi = Math.min(z2, L.bp);
      if (hi <= lo) continue;
      if (lo > z) M = mul(FL.drift(lo - z), M);
      const k2 = k2Of(id, powers[id]);
      M = mul(k2 > 0 ? L.matrix(k2, lo, hi) : FL.drift(hi - lo), M);
      z = hi;
    }
    if (z2 > z) M = mul(FL.drift(z2 - z), M);
    return M;
  }
  /** Larmor rotation, radians, accumulated between z1 and z2 (all lenses taken with the same polarity). */
  function rotation(powers, z1, z2) {
    let r = 0;
    for (const id of order) {
      const L = lens[id], lo = Math.max(z1, L.bm), hi = Math.min(z2, L.bp);
      if (hi > lo) r += L.rotation(k2Of(id, powers[id]), lo, hi);
    }
    return r;
  }

  /** A sampled trajectory from (z0, y0, u0) to zEnd: arrays z, y, u (the points to draw, dense inside each
   *  field and just the ends of a drift), plus an exact at(z) and the axis crossings. `perLens` samples per
   *  field. The fan of rays on the page is a handful of multiples of one of these. */
  function trace(powers, z0, y0, u0, zEnd, perLens) {
    const n = perLens || 48;
    const zs = [z0], ys = [y0], us = [u0];
    const segs = [];                                      // for refining crossings inside a field
    let z = z0, y = y0, u = u0;
    const step = (M, ze) => { const ny = M.A * y + M.B * u, nu = M.C * y + M.D * u; y = ny; u = nu; z = ze; };
    const push = () => { zs.push(z); ys.push(y); us.push(u); };
    for (const id of order) {
      const L = lens[id], lo = Math.max(z, L.bm), hi = Math.min(zEnd, L.bp);
      if (hi <= lo) continue;
      if (lo > z) { step(FL.drift(lo - z), lo); push(); }
      const k2 = k2Of(id, powers[id]);
      if (k2 > 0) {
        const g = FL.grid(L, lo, hi, n), y0s = y, u0s = u;
        segs.push({ kind: 'field', L: L, k2: k2, z1: lo, z2: hi, y0: y0s, u0: u0s });
        for (let j = 1; j < g.length; j++) {
          const M = L.matrix(k2, lo, g[j]);
          zs.push(g[j]); ys.push(M.A * y0s + M.B * u0s); us.push(M.C * y0s + M.D * u0s);
        }
        const Mf = L.matrix(k2, lo, hi); y = Mf.A * y0s + Mf.B * u0s; u = Mf.C * y0s + Mf.D * u0s; z = hi;
      } else {
        step(FL.drift(hi - lo), hi); push();
      }
    }
    if (zEnd > z) { step(FL.drift(zEnd - z), zEnd); push(); }

    const at = zq => { const M = chain(powers, z0, zq); return { y: M.A * y0 + M.B * u0, u: M.C * y0 + M.D * u0 }; };
    const zeros = [];
    for (let i = 0; i + 1 < zs.length; i++) {
      const ya = ys[i], yb = ys[i + 1];
      if (ya === 0 || !(ya * yb <= 0) || zs[i + 1] === zs[i]) continue;
      // bracket [zs[i], zs[i+1]]: straight inside a drift, refined by bisection on the exact at() elsewhere
      let lo = zs[i], hi = zs[i + 1], flo = ya;
      if (us[i] !== 0 && Math.abs(yb - (ya + us[i] * (hi - lo))) < 1e-12 * (Math.abs(ya) + Math.abs(yb) + 1e-30)) {
        const zc = lo - ya / us[i];
        if (zc > z0 + 1e-9) zeros.push(zc);
        continue;
      }
      for (let k = 0; k < 60; k++) {
        const mid = (lo + hi) / 2, ym = at(mid).y;
        if (flo * ym <= 0) hi = mid; else { lo = mid; flo = ym; }
      }
      const zc = (lo + hi) / 2;
      if (zc > z0 + 1e-9) zeros.push(zc);
    }
    return { z: zs, y: ys, u: us, at: at, zeros: zeros, yEnd: y, uEnd: u };
  }

  /* ------------------------------------------------------------ the effective thin-lens column */
  /** Replace every thick lens by the thin lens that is exactly equivalent to it, for `powers`.
   *
   *  Any 2x2 ray matrix with C != 0 is drift(dOut) . thin(f) . drift(dIn) (FieldLens.cardinal). Doing that for
   *  each lens piece -- the objective in two, cut at the reference plane -- gives a column of ideal thin lenses
   *  that propagates rays between the pieces' boundary planes EXACTLY. dIn + dOut is not the piece's length
   *  (a thick lens can shorten or lengthen the column), so effective positions are the real ones minus the
   *  running total of those differences, with the reference plane anchored at its own z.
   *
   *  Returns { lenses: [{id, z, power}] (OL appears as OL1 and OL2), toEff(z), planes } */
  function effective(powers, zRef) {
    const pieces = [];
    for (const id of order) {
      const L = lens[id], k2 = k2Of(id, powers[id]);
      const cut = L.kind === 'half' ? [[L.bm, zRef, 'OL1'], [zRef, L.bp, 'OL2']] : [[L.bm, L.bp, id]];
      for (const c of cut) {
        const len = c[1] - c[0];
        const card = k2 > 0 ? FL.cardinal(L.matrix(k2, c[0], c[1]), len) : { f: Infinity, dIn: len / 2, dOut: len / 2 };
        pieces.push({ id: c[2], z1: c[0], z2: c[1], e: c[0] + card.dIn, delta: len - card.dIn - card.dOut,
                      power: isFinite(card.f) ? 1 / card.f : 0 });
      }
    }
    // running total of length differences strictly upstream of each piece / plane
    const before = z => pieces.reduce((s, p) => s + (p.z2 <= z + 1e-9 ? p.delta : 0), 0);
    const c0 = before(zRef);                              // anchor: the reference plane keeps its own z
    const toEff = z => z - before(z) + c0;
    const lensesEff = pieces.map(p => {
      const upstream = pieces.reduce((s, q) => s + (q.z2 <= p.z1 + 1e-9 ? q.delta : 0), 0);
      return { id: p.id, z: p.e - upstream + c0, power: p.power };
    });
    return { lenses: lensesEff, toEff: toEff,
             planes: { SRC: toEff(Z_SRC), VOA: toEff(Z_VOA), REF: zRef, EAP: toEff(Z_EAP), EFP: toEff(Z_EFP) } };
  }

  return { Z_MIN, Z_MAX, Z_SRC, Z_VOA, Z_OL, Z_REF, Z_EAP, Z_EFP, VOA_R, LENSES, F_MIN, CAP, WIDTH, lens, order,
           k2Of, pMax, chain, rotation, trace, effective };
})();

if (typeof window !== 'undefined') window.ElectronColumn = ElectronColumn;
