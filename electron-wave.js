/* ================================================================
   electron-wave.js -- the electron as a WAVE on the column's paraxial optics.

   Pure logic plus a painter; no page state. The page hands it the exact ray-transfer matrix from the specimen
   to a plane (ElectronColumn.chain, through the real lens fields) and the Larmor rotation accumulated on the
   way, and gets back the intensity there.

   The probe is the ideal converging wave the ray model describes, with the aberration phase in Krivanek notation
   up to third order; the specimen is a thin phase object; the plane is reached with the Collins integral, the
   wave-optical counterpart of the ABCD matrix. Everything is in the rotating frame, like the rays, and the image
   is turned back by the Larmor angle when it is drawn.
   ================================================================ */
"use strict";

const WaveImages = (function () {

  /* ------------------------------------------------------------ physics constants */
  // CODATA 2018 exact / recommended values
  const H_PLANCK = 6.62607015e-34, M_E = 9.1093837015e-31, Q_E = 1.602176634e-19, C_LIGHT = 299792458;
  /** Relativistic electron wavelength in nm for an accelerating voltage in kV. */
  function wavelengthNm(kV) {
    const V = kV * 1e3;
    const p = Math.sqrt(2 * M_E * Q_E * V * (1 + Q_E * V / (2 * M_E * C_LIGHT * C_LIGHT)));
    return H_PLANCK / p * 1e9;
  }

  /** Aberration coefficients, Krivanek notation, nm. C10, C30 are round; the rest have a (cos) and b (sin)
   *  components. The Nion residuals are the order of magnitude a tuned corrector leaves: ~2 nm at first order,
   *  ~20 nm at second, ~200 nm at third. */
  const ABERR = [
    { k: 'C10', n: 1, m: 0 }, { k: 'C12', n: 1, m: 2 },
    { k: 'C21', n: 2, m: 1 }, { k: 'C23', n: 2, m: 3 },
    { k: 'C30', n: 3, m: 0 }, { k: 'C32', n: 3, m: 2 }, { k: 'C34', n: 3, m: 4 },
  ];
  // a tuned corrector leaves ~2 nm first order, ~20 nm second, ~200 nm third, in both components (a and b) of each
  // non-axial term; C10 and C30 are axial and have no b
  const NION_RESIDUALS = { C10: [0, 0], C12: [2, 2], C21: [20, 20], C23: [20, 20],
                           C30: [200, 0], C32: [200, 200], C34: [200, 200] };

  /* ------------------------------------------------------------ FFT: iterative radix-2, in place */
  const _plans = {};
  function plan(n) {
    if (_plans[n]) return _plans[n];
    const bits = Math.round(Math.log2(n));
    const rev = new Uint32Array(n);
    for (let i = 0; i < n; i++) {
      let r = 0;
      for (let b = 0; b < bits; b++) r |= ((i >> b) & 1) << (bits - 1 - b);
      rev[i] = r;
    }
    const cos = new Float64Array(n / 2), sin = new Float64Array(n / 2);
    for (let k = 0; k < n / 2; k++) { cos[k] = Math.cos(2 * Math.PI * k / n); sin[k] = Math.sin(2 * Math.PI * k / n); }
    return (_plans[n] = { rev, cos, sin });
  }
  /** Unnormalised DFT, exp(-2 pi i jk/n) forward and exp(+) inverse. */
  function fft1(xr, xi, n, inv) {
    const P = plan(n), rev = P.rev, sg = inv ? 1 : -1;
    for (let i = 0; i < n; i++) {
      const j = rev[i];
      if (j > i) { let t = xr[i]; xr[i] = xr[j]; xr[j] = t; t = xi[i]; xi[i] = xi[j]; xi[j] = t; }
    }
    for (let size = 2; size <= n; size <<= 1) {
      const half = size >> 1, step = n / size;
      for (let i = 0; i < n; i += size)
        for (let k = 0; k < half; k++) {
          const wr = P.cos[k * step], wi = sg * P.sin[k * step];
          const a = i + k, b = a + half;
          const tr = xr[b] * wr - xi[b] * wi, ti = xr[b] * wi + xi[b] * wr;
          xr[b] = xr[a] - tr; xi[b] = xi[a] - ti;
          xr[a] += tr;        xi[a] += ti;
        }
    }
  }
  function fft2(re, im, N, inv) {
    for (let r = 0; r < N; r++) fft1(re.subarray(r * N, r * N + N), im.subarray(r * N, r * N + N), N, inv);
    const cr = new Float64Array(N), ci = new Float64Array(N);
    for (let c = 0; c < N; c++) {
      for (let r = 0; r < N; r++) { cr[r] = re[r * N + c]; ci[r] = im[r * N + c]; }
      fft1(cr, ci, N, inv);
      for (let r = 0; r < N; r++) { re[r * N + c] = cr[r]; im[r * N + c] = ci[r]; }
    }
  }
  /** Swap quadrants, so index N/2 is the origin. Its own inverse for even N. */
  function shift2(re, im, N) {
    const h = N / 2;
    for (let i = 0; i < h; i++)
      for (let j = 0; j < N; j++) {
        const a = i * N + j, b = (i + h) * N + ((j + h) & (N - 1));
        let t = re[a]; re[a] = re[b]; re[b] = t;
        t = im[a]; im[a] = im[b]; im[b] = t;
      }
  }
  /** Centred FFT: the arrays hold the origin at index N/2 on the way in and on the way out. */
  function cfft2(re, im, N, inv) { shift2(re, im, N); fft2(re, im, N, inv); shift2(re, im, N); }

  /* ------------------------------------------------------------ the simulation.
     Units: specimen plane in nm, angles in rad, column lengths (ABCD) in mm converted to nm where they meet lambda. */

  /** Pixel size and grid for this render. The pixel must resolve every angle in the wave -- the probe cone plus two
   *  Bragg orders for a lattice, or the scattering of the amorphous grain -- times the oversample factor S. For a
   *  lattice the window is then nudged to a whole number of periods, so the FFT's periodic boundary sees a
   *  seamless crystal. */
  function chooseGrid(p) {
    const lam = p.lambda;
    let extra = 0;
    if (p.spec.type === 'lattice') extra = 2 * lam / (p.spec.dA / 10);
    else if (p.spec.type === 'amorph') extra = 1.5 * lam / p.spec.ellNm;
    const thetaMax = Math.max(p.S * (p.alpha + extra), p.collimated ? 1e-4 : p.S * p.alpha, 1e-5);
    let dx = lam / (2 * thetaMax);
    if (p.spec.type === 'lattice') {
      const d = p.spec.dA / 10, m = Math.max(1, Math.round(p.N * dx / d));
      dx = m * d / p.N;
    }
    return { N: p.N, dx: dx };
  }

  function chi(p, tx, ty) {
    const th = Math.sqrt(tx * tx + ty * ty), ph = Math.atan2(ty, tx);
    const a = p.aberr;
    let s = 0;
    for (const c of ABERR) {
      const [ca, cb] = a[c.k] || [0, 0];
      const v = c.k === 'C10' ? ca + p.defocusNm : ca;
      const ang = c.m === 0 ? v : v * Math.cos(c.m * ph) + cb * Math.sin(c.m * ph);
      if (!ang) continue;
      s += ang * Math.pow(th, c.n + 1) / (c.n + 1);
    }
    return 2 * Math.PI / p.lambda * s;
  }

  function normalise(re, im) {
    let s = 0;
    for (let k = 0; k < re.length; k++) s += re[k] * re[k] + im[k] * im[k];
    const f = s > 0 ? 1 / Math.sqrt(s) : 0;
    for (let k = 0; k < re.length; k++) { re[k] *= f; im[k] *= f; }
  }

  /** Incident wave at the specimen, centred, normalised to unit total intensity. */
  function buildProbe(p, g) {
    const N = g.N, dx = g.dx, NN = N * N;
    const re = new Float64Array(NN), im = new Float64Array(NN);
    if (p.collimated) {                     // plane wave; the illuminated disc is far wider than the window
      re.fill(1 / N);
      return { re, im };
    }
    const dk = 1 / (N * dx), dth = p.lambda * dk;
    for (let i = 0; i < N; i++) {
      const ty = (i - N / 2) * dth;
      for (let j = 0; j < N; j++) {
        const tx = (j - N / 2) * dth, th = Math.hypot(tx, ty);
        // one-pixel soft edge, so the disc's area is right to a fraction of a pixel
        const amp = Math.min(1, Math.max(0, (p.alpha - th) / dth + 0.5));
        if (amp <= 0) continue;
        const c = chi(p, tx, ty);
        re[i * N + j] = amp * Math.cos(c); im[i * N + j] = -amp * Math.sin(c);
      }
    }
    cfft2(re, im, N, true);
    normalise(re, im);
    return { re, im };
  }

  /** Seeded PRNG, so a given amorphous specimen is the same film from render to render. */
  function mulberry32(a) {
    return () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a);
                   t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
  }
  /** Specimen phase phi(x), rad, or null for vacuum. */
  function specimenPhase(p, g) {
    const N = g.N, dx = g.dx, NN = N * N;
    if (p.spec.type === 'vacuum') return null;
    const ph = new Float64Array(NN);
    if (p.spec.type === 'lattice') {
      const d = p.spec.dA / 10, s = p.spec.sigma;
      for (let i = 0; i < N; i++) {
        const cy = Math.cos(2 * Math.PI * (i - N / 2) * dx / d);
        for (let j = 0; j < N; j++) ph[i * N + j] = s * (Math.cos(2 * Math.PI * (j - N / 2) * dx / d) + cy);
      }
      return ph;
    }
    // amorphous: white noise, Gaussian-filtered to the grain size, scaled to the requested rms
    const rnd = mulberry32(12345), re = new Float64Array(NN), im = new Float64Array(NN);
    for (let k = 0; k < NN; k++) {                        // Box-Muller
      re[k] = Math.sqrt(-2 * Math.log(rnd() + 1e-300)) * Math.cos(2 * Math.PI * rnd());
    }
    cfft2(re, im, N, false);
    const dk = 1 / (N * dx), l = p.spec.ellNm;
    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++) {
        const k2 = ((i - N / 2) ** 2 + (j - N / 2) ** 2) * dk * dk;
        const f = Math.exp(-Math.PI * Math.PI * k2 * l * l);
        re[i * N + j] *= f; im[i * N + j] *= f;
      }
    cfft2(re, im, N, true);
    let s = 0;
    for (let k = 0; k < NN; k++) s += re[k] * re[k];
    const f = p.spec.sigma / Math.sqrt(s / NN || 1);
    for (let k = 0; k < NN; k++) ph[k] = re[k] * f;
    return ph;
  }

  /** Propagate the exit wave through the paraxial transfer M = {A, B} (B in mm) to the plane.
   *
   *  Collins integral, psi2(u) ~ exp(i pi D u^2/(lam B)) * FT[ psi1(x) exp(i pi A x^2/(lam B)) ](u/(lam B)).
   *  Two evaluations, each exact where it samples well:
   *    |B/A| >= z*   chirp, then one FFT; output pixel lam|B|/(N dx)     (diffraction-like planes)
   *    |B/A| <  z*   angular-spectrum propagation over B/A, then magnify by A; pixel |A| dx
   *                  (the decomposition M = [[A,0],[C,1/A]] . [[1,B/A],[0,1]])     (image-like planes)
   *  with z* = N dx^2 / lam, which is exactly where the first one's pre-chirp and the second one's transfer
   *  function each stop being sampled at Nyquist. The outgoing phase (the D term, and the 1/A lens) never changes an
   *  intensity, so it is not computed.
   *  Returns { I, du (mm), branch, flip } -- flip is true when the plane is inverted (A or B < 0). */
  function propagateField(re0, im0, g, lam, M) {
    const N = g.N, dx = g.dx;
    const re = Float64Array.from(re0), im = Float64Array.from(im0);
    const A = M.A, Bnm = M.B * 1e6, zStar = N * dx * dx / lam;
    let delta, branch, flip;
    if (Math.abs(Bnm) >= Math.abs(A) * zStar) {
      branch = 'fourier';
      const q = Math.PI * A / (lam * Bnm);
      for (let i = 0; i < N; i++) {
        const y = (i - N / 2) * dx;
        for (let j = 0; j < N; j++) {
          const x = (j - N / 2) * dx, ph = q * (x * x + y * y), k = i * N + j;
          const c = Math.cos(ph), s = Math.sin(ph), r = re[k], m = im[k];
          re[k] = r * c - m * s; im[k] = r * s + m * c;
        }
      }
      cfft2(re, im, N, false);
      delta = lam * Bnm / (N * dx);                     // signed nm per output pixel: u = (j - N/2) delta
      flip = Bnm < 0;
    } else {
      branch = 'fresnel';
      const z = Bnm / A, dk = 1 / (N * dx);
      cfft2(re, im, N, false);
      for (let i = 0; i < N; i++) {
        const ky = (i - N / 2) * dk;
        for (let j = 0; j < N; j++) {
          const kx = (j - N / 2) * dk, ph = -Math.PI * lam * z * (kx * kx + ky * ky), k = i * N + j;
          const c = Math.cos(ph), s = Math.sin(ph), r = re[k], m = im[k];
          re[k] = r * c - m * s; im[k] = r * s + m * c;
        }
      }
      cfft2(re, im, N, true);
      delta = A * dx;                                   // u = (j - N/2) delta (the image is magnified by A)
      flip = A < 0;
    }
    return { re, im, delta, branch, flip, zStar };
  }
  function propagate(re0, im0, g, lam, M) {
    const f = propagateField(re0, im0, g, lam, M), NN = g.N * g.N;
    const I = new Float64Array(NN);
    let s = 0;
    for (let k = 0; k < NN; k++) { I[k] = f.re[k] * f.re[k] + f.im[k] * f.im[k]; s += I[k]; }
    if (s > 0) for (let k = 0; k < NN; k++) I[k] /= s;
    return { I, du: Math.abs(f.delta) * 1e-6, branch: f.branch, flip: f.flip, zStar: f.zStar };
  }

  /** Full width at half maximum along x through the brightest pixel, linear interpolation, nm. */
  function fwhmOf(I, g) {
    const N = g.N;
    let kMax = 0;
    for (let k = 1; k < I.length; k++) if (I[k] > I[kMax]) kMax = k;
    const r = Math.floor(kMax / N), c = kMax % N, half = I[kMax] / 2;
    const at = j => I[r * N + j];
    let j1 = c, j2 = c;
    while (j1 > 0 && at(j1 - 1) >= half) j1--;
    while (j2 < N - 1 && at(j2 + 1) >= half) j2++;
    if (j1 === 0 || j2 === N - 1) return NaN;
    const left  = (j1 - 1) + (half - at(j1 - 1)) / (at(j1) - at(j1 - 1));
    const right = j2 + (at(j2) - half) / (at(j2) - at(j2 + 1));
    return (right - left) * g.dx;
  }

  /* ------------------------------------------------------------ finite specimen height (dz != 0)

     A probe focused on the crossover plane z_c is a DISC of radius dz*alpha on a specimen dz away. Its wave-optical
     description needs a grid that holds that disc AND resolves its defocus phase (the local angle reaches alpha), which
     grows with dz and is out of reach beyond a micron or so. The way out is to never put the wave on the specimen.

     A thin periodic specimen is a sum of plane-wave modulations, T(r) = sum_g t_g exp(2 pi i g.r). Through any ABCD
     system with matrix Ms (specimen -> plane), modulation by exp(2 pi i g.r) only SHIFTS the output and adds a phase:

         out_g(u) = exp[2 pi i Ds g.u - i pi Ds Bs lam g^2] * F(u - lam Bs g),

     where F is the output for the UNMODULATED wave. F is the same function whatever the specimen, and by the composition
     of ABCD systems it equals the propagation of the focused probe from its crossover through Mc = Ms . Mdz, which is
     the well-sampled dz = 0 calculation. So the pattern is   out(u) = sum_g t_g out_g(u)   -- a handful of copies of F,
     displaced by lam Bs g, with the exact relative phases. Nothing in this depends on how far the specimen is.

     What each plane shows follows from F. At a diffraction-like plane F is the pupil disc (radius |Bc| alpha), so the
     beams are discs displaced by lam Bs g and they interfere where they overlap. At an image-like plane F is the
     magnified focused probe, a spot, so the beams are separate spots displaced by lam Bs g = A dz theta_g -- the
     image-coupled aperture shows SPOTS however far the specimen is from the crossover. */

  /** Bessel function J_n(x) of integer order by its power series (fine for |x| up to ~20). */
  function besselJ(n, x) {
    const an = Math.abs(n), h = x / 2;
    let term = an === 0 ? 1 : Math.pow(h, an);
    for (let k = 2; k <= an; k++) term /= k;
    let s = term;
    for (let k = 1; k < 120; k++) {
      term *= -h * h / (k * (k + an));
      s += term;
      if (k > 3 && Math.abs(term) < 1e-18 * Math.abs(s)) break;
    }
    return (n < 0 && an % 2) ? -s : s;
  }

  /** Bragg beams of the lattice phase object exp[i s (cos(2 pi x/d) + cos(2 pi y/d))]:
   *  exp(i s cos u) = sum_n i^n J_n(s) exp(i n u), so t(m,n) = i^(m+n) J_m(s) J_n(s) at g = (m, n)/d. */
  function latticeBeams(p) {
    const d = p.spec.dA / 10, s = p.spec.sigma, J = [];
    let nmax = 0;
    for (let n = 0; n <= 16; n++) { J[n] = besselJ(n, s); if (Math.abs(J[n]) > 1e-4) nmax = n; }
    const jn = n => (n < 0 && (-n) % 2) ? -J[-n] : J[Math.abs(n)];
    const PH = [[1, 0], [0, 1], [-1, 0], [0, -1]], out = [];
    for (let m = -nmax; m <= nmax; m++)
      for (let n = -nmax; n <= nmax; n++) {
        const w = jn(m) * jn(n);
        if (w * w < 1e-5 && !(m === 0 && n === 0)) continue;     // below the display's four decades
        const ph = PH[(((m + n) % 4) + 4) % 4];
        out.push({ gx: m / d, gy: n / d, re: w * ph[0], im: w * ph[1] });
      }
    return out;
  }

  /** The strongest `cap` Fourier components of exp(i phi) for a periodic phase screen on the grid (amorphous film). */
  function screenBeams(p, g, cap) {
    const N = g.N, NN = N * N, phase = specimenPhase(p, g);
    if (!phase) return [{ gx: 0, gy: 0, re: 1, im: 0 }];
    const re = new Float64Array(NN), im = new Float64Array(NN);
    for (let k = 0; k < NN; k++) { re[k] = Math.cos(phase[k]); im[k] = Math.sin(phase[k]); }
    cfft2(re, im, N, false);
    const idx = [];
    for (let k = 0; k < NN; k++) { const w = (re[k] * re[k] + im[k] * im[k]) / (NN * NN); if (w > 1e-9) idx.push([w, k]); }
    idx.sort((a, b) => b[0] - a[0]);
    const dk = 1 / (N * g.dx), out = [];
    for (const [, k] of idx.slice(0, cap)) {
      const i = Math.floor(k / N), j = k % N;
      out.push({ gx: (j - N / 2) * dk, gy: (i - N / 2) * dk, re: re[k] / NN, im: im[k] / NN });
    }
    return out;
  }

  const sinc = x => Math.abs(x) < 1e-9 ? 1 : Math.sin(x) / x;

  /** The dz != 0 render: see the block comment above. Returns the same shape as simulate(). p.Mc is the matrix from
   *  the probe's crossover to the plane; p.M is the one from the specimen. */
  function simulateDefocused(p) {
    const lam = p.lambda, p0 = Object.assign({}, p, { defocusNm: 0 });
    const g = chooseGrid(p0), N = g.N, NN = N * N, win = N * g.dx;
    const probe = buildProbe(p0, g);
    const I1 = new Float64Array(NN);
    for (let k = 0; k < NN; k++) I1[k] = probe.re[k] * probe.re[k] + probe.im[k] * probe.im[k];
    const F = propagateField(probe.re, probe.im, g, lam, p.Mc);
    const Bs = p.M.B * 1e6, Ds = p.M.D, Bc = p.Mc.B * 1e6, Dc = p.Mc.D, Ac = p.Mc.A, Cc = p.Mc.C * 1e-6;
    const coherentOK = p.spec.type !== 'amorph';
    const beams = p.spec.type === 'lattice' ? latticeBeams(p)
                : p.spec.type === 'amorph'  ? screenBeams(p, g, 250)
                : [{ gx: 0, gy: 0, re: 1, im: 0 }];
    let pw = 0; for (const b of beams) pw += b.re * b.re + b.im * b.im;
    for (const b of beams) { b.t2 = b.re * b.re + b.im * b.im; b.sx = lam * Bs * b.gx; b.sy = lam * Bs * b.gy; }
    const Iout = new Float64Array(NN), delta = F.delta;
    let branch = F.branch, du = Math.abs(delta) * 1e-6, flip = F.flip;

    if (F.branch === 'fourier') {
      // ---- diffraction-like plane: discs displaced by lam Bs g, added with their exact phases
      const qc = Math.PI * Dc / (lam * Bc);
      // support of F, to skip the empty plane quickly
      let i0 = N, i1 = -1, j0 = N, j1 = -1, mx = 0;
      const PhiI = new Float64Array(NN);
      for (let k = 0; k < NN; k++) { PhiI[k] = F.re[k] * F.re[k] + F.im[k] * F.im[k]; if (PhiI[k] > mx) mx = PhiI[k]; }
      for (let i = 0; i < N; i++) for (let j = 0; j < N; j++)
        if (PhiI[i * N + j] > 1e-6 * mx) { if (i < i0) i0 = i; if (i > i1) i1 = i; if (j < j0) j0 = j; if (j > j1) j1 = j; }
      const nb = beams.length;
      for (const b of beams) {
        b.ox = b.sx / delta; b.oy = b.sy / delta;                               // pixel offsets of the copies
        b.bi0 = Math.floor(i0 + b.oy) - 1; b.bi1 = Math.ceil(i1 + b.oy) + 1;       // where this copy lives
        b.bj0 = Math.floor(j0 + b.ox) - 1; b.bj1 = Math.ceil(j1 + b.ox) + 1;
      }
      // pair factors: pixel-averaged cross terms (the phase difference of two beams is linear in u)
      const S = new Float64Array(nb * nb);
      if (coherentOK)
        for (let a = 0; a < nb; a++) for (let c = 0; c < nb; c++) {
          const Kx = 2 * Math.PI * Ds * (beams[a].gx - beams[c].gx) - 2 * qc * (beams[a].sx - beams[c].sx);
          const Ky = 2 * Math.PI * Ds * (beams[a].gy - beams[c].gy) - 2 * qc * (beams[a].sy - beams[c].sy);
          S[a * nb + c] = sinc(Kx * Math.abs(delta) / 2) * sinc(Ky * Math.abs(delta) / 2);
        }
      const cr = new Float64Array(nb), ci = new Float64Array(nb), act = new Int32Array(nb);
      const PHr = F.re, PHi = F.im;
      let ubi0 = N, ubi1 = -1, ubj0 = N, ubj1 = -1;
      for (const b of beams) {
        ubi0 = Math.min(ubi0, Math.floor(i0 + b.oy) - 1); ubi1 = Math.max(ubi1, Math.ceil(i1 + b.oy) + 1);
        ubj0 = Math.min(ubj0, Math.floor(j0 + b.ox) - 1); ubj1 = Math.max(ubj1, Math.ceil(j1 + b.ox) + 1);
      }
      ubi0 = Math.max(0, ubi0); ubj0 = Math.max(0, ubj0); ubi1 = Math.min(N - 1, ubi1); ubj1 = Math.min(N - 1, ubj1);
      for (let i = ubi0; i <= ubi1; i++) {
        const uy = (i - N / 2) * delta;
        for (let j = ubj0; j <= ubj1; j++) {
          const ux = (j - N / 2) * delta;
          let na = 0;
          for (let q = 0; q < nb; q++) {
            const b = beams[q];
            if (i < b.bi0 || i > b.bi1 || j < b.bj0 || j > b.bj1) continue;
            const fi = i - b.oy, fj = j - b.ox;
            if (fi < 0 || fj < 0 || fi > N - 1 || fj > N - 1) continue;
            const a0 = Math.floor(fi), b0 = Math.floor(fj), ty = fi - a0, tx = fj - b0;
            const a1 = Math.min(N - 1, a0 + 1), b1 = Math.min(N - 1, b0 + 1);
            const k00 = a0 * N + b0, k01 = a0 * N + b1, k10 = a1 * N + b0, k11 = a1 * N + b1;
            const w00 = (1 - ty) * (1 - tx), w01 = (1 - ty) * tx, w10 = ty * (1 - tx), w11 = ty * tx;
            const fr = w00 * PHr[k00] + w01 * PHr[k01] + w10 * PHr[k10] + w11 * PHr[k11];
            const fm = w00 * PHi[k00] + w01 * PHi[k01] + w10 * PHi[k10] + w11 * PHi[k11];
            if (fr === 0 && fm === 0) continue;
            const dxu = ux - b.sx, dyu = uy - b.sy;
            const th = 2 * Math.PI * Ds * (b.gx * ux + b.gy * uy) - Math.PI * Ds * Bs * lam * (b.gx * b.gx + b.gy * b.gy)
                     + qc * (dxu * dxu + dyu * dyu);
            const cs = Math.cos(th), sn = Math.sin(th);
            const ar = fr * cs - fm * sn, ai = fr * sn + fm * cs;           // F * exp(i theta)
            cr[na] = b.re * ar - b.im * ai; ci[na] = b.re * ai + b.im * ar; act[na] = q; na++;
          }
          if (!na) continue;
          let v = 0;
          for (let a = 0; a < na; a++) v += cr[a] * cr[a] + ci[a] * ci[a];
          if (coherentOK)
            for (let a = 0; a < na; a++)
              for (let c = a + 1; c < na; c++)
                v += 2 * S[act[a] * nb + act[c]] * (cr[a] * cr[c] + ci[a] * ci[c]);
          Iout[i * N + j] = v;
        }
      }
    } else {
      // ---- image-like plane: the beams are well-separated magnified copies of the focused probe, i.e. spots. Their
      // separations (lam Bs g) dwarf a spot, so they add incoherently; draw each as a Gaussian of the probe's width on
      // a window that holds them all.
      branch = 'spots';
      const PhiI = new Float64Array(NN);
      for (let k = 0; k < NN; k++) PhiI[k] = F.re[k] * F.re[k] + F.im[k] * F.im[k];
      let fw = fwhmOf(PhiI, g) * Math.abs(Ac);
      if (!isFinite(fw) || fw <= 0) fw = 2 * Math.abs(delta);
      const sig = fw / 2.3548;
      let reach = 0;
      for (const b of beams) if (b.t2 > 1e-5) reach = Math.max(reach, Math.hypot(b.sx, b.sy));
      const W = Math.max(2.4 * reach, 40 * sig, 16 * sig);
      const Dlt = W / N, se2 = sig * sig + 1.7 * Dlt * Dlt,        // drawn about 1.3 display pixels wide so they can be seen;
        // the flux of each spot is exact, only its picture is larger than the 1.4 nm probe
        se = Math.sqrt(se2), sA = Ac < 0 ? -1 : 1;
      const rp = Math.ceil(4 * se / Dlt);
      for (const b of beams) {
        const cx = N / 2 + b.sx * sA / Dlt, cy = N / 2 + b.sy * sA / Dlt;
        if (cx < -rp || cy < -rp || cx > N + rp || cy > N + rp) continue;
        const amp = b.t2 / pw / (2 * Math.PI * se2) * Dlt * Dlt;
        for (let i = Math.max(0, Math.round(cy) - rp); i <= Math.min(N - 1, Math.round(cy) + rp); i++)
          for (let j = Math.max(0, Math.round(cx) - rp); j <= Math.min(N - 1, Math.round(cx) + rp); j++) {
            const ex = (j - cx) * Dlt, ey = (i - cy) * Dlt;
            Iout[i * N + j] += amp * Math.exp(-(ex * ex + ey * ey) / (2 * se2));
          }
      }
      du = Dlt * 1e-6;
    }
    let s = 0;
    for (let k = 0; k < NN; k++) s += Iout[k];
    if (s > 0) for (let k = 0; k < NN; k++) Iout[k] /= s;
    const out = { I: Iout, du: du, branch: branch, flip: flip, zStar: F.zStar };
    return { g, win, I1, phase: null, M: p.M, rot: p.rot || 0, out, outOfWindow: false, probeR: 0, defocused: true,
             nBeams: beams.length, coherent: coherentOK && F.branch === 'fourier', fwhm: fwhmOf(I1, g) };
  }

  /** One full render's worth of physics. p carries everything, including p.M (the exact transfer matrix from the
   *  specimen to the plane) and p.rot (the Larmor rotation on the way); nothing here reads the page. */
  function simulate(p) {
    const g = chooseGrid(p);
    const win = g.N * g.dx;                                  // nm
    // Is the probe inside the window? Geometric radius |dz| alpha, plus a few diffraction widths.
    const probeR = p.collimated ? 0 : Math.abs(p.defocusNm + (p.aberr.C10 ? p.aberr.C10[0] : 0)) * p.alpha
                                    + 3 * p.lambda / Math.max(p.alpha, 1e-9);
    const outOfWindow = !p.collimated && probeR > 0.45 * win;
    if (!p.collimated && p.Mc && (outOfWindow || p.forceBragg)) return simulateDefocused(p);
    const probe = buildProbe(p, g);
    const phase = specimenPhase(p, g);
    const exRe = Float64Array.from(probe.re), exIm = Float64Array.from(probe.im);
    if (phase)
      for (let k = 0; k < exRe.length; k++) {
        const c = Math.cos(phase[k]), s = Math.sin(phase[k]), r = exRe[k], m = exIm[k];
        exRe[k] = r * c - m * s; exIm[k] = r * s + m * c;
      }
    const I1 = new Float64Array(g.N * g.N);
    for (let k = 0; k < I1.length; k++) I1[k] = probe.re[k] * probe.re[k] + probe.im[k] * probe.im[k];
    const out = propagate(exRe, exIm, g, p.lambda, p.M);
    return { g, win, I1, phase, M: p.M, rot: p.rot || 0, out, outOfWindow, probeR,
             fwhm: p.collimated ? NaN : fwhmOf(I1, g) };
  }

  /* ------------------------------------------------------------ drawing */
  function lut(name, t) {
    t = Math.min(1, Math.max(0, t));
    if (name === 'gray') { const v = Math.round(255 * t); return [v, v, v]; }
    return [Math.round(255 * Math.min(1, 3 * t)), Math.round(255 * Math.min(1, Math.max(0, 3 * t - 1))),
            Math.round(255 * Math.min(1, Math.max(0, 3 * t - 2)))];
  }
  /** Paint an N x N intensity onto a canvas. opt: log/dyn (log scale over `dyn` decades), signed, cmap, flip (an
   *  inverted plane is mirrored back upright), rot (radians; the picture is turned by this about its centre). */
  function paint(cv, I, N, opt) {
    const M = opt.size || N;                    // canvas pixels: the data are N x N, stretched smoothly to M x M
    cv.width = M; cv.height = M;
    const g = cv.getContext('2d'), img = g.createImageData(N, N);
    let mx = 0, mn = Infinity;
    for (let k = 0; k < I.length; k++) { if (I[k] > mx) mx = I[k]; if (I[k] < mn) mn = I[k]; }
    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++) {
        // canvas row 0 is the top: draw +y up, and mirror an inverted plane back into place
        const si = opt.flip ? i : N - 1 - i, sj = opt.flip ? N - 1 - j : j;
        const v = I[si * N + sj];
        let t;
        if (opt.signed) t = mx > mn ? (v - mn) / (mx - mn) : 0.5;
        else if (opt.log) t = v > 0 && mx > 0 ? 1 + Math.log10(v / mx) / opt.dyn : 0;
        else t = mx > 0 ? v / mx : 0;
        if (opt.light) t = 1 - t;                 // light theme: zero intensity is white
        const c = lut(opt.cmap, t), o = (i * N + j) * 4;
        img.data[o] = c[0]; img.data[o + 1] = c[1]; img.data[o + 2] = c[2]; img.data[o + 3] = 255;
      }
    if (opt.rot || M !== N) {
      const tmp = document.createElement('canvas'); tmp.width = N; tmp.height = N;
      tmp.getContext('2d').putImageData(img, 0, 0);
      g.fillStyle = opt.light ? '#fff' : '#000'; g.fillRect(0, 0, M, M);
      if (opt.rot) {
        // the turned square leaves corners nobody computed: hatch them so they cannot be mistaken for zero intensity
        g.fillStyle = opt.light ? '#e6eaed' : '#16161b'; g.fillRect(0, 0, M, M);
        g.strokeStyle = opt.light ? '#d0d5da' : '#2e2e38'; g.lineWidth = Math.max(1, M / 256);
        for (let d = -M; d < M; d += Math.max(6, M / 42)) { g.beginPath(); g.moveTo(d, 0); g.lineTo(d + M, M); g.stroke(); }
      }
      g.save(); g.translate(M / 2, M / 2); g.rotate(-(opt.rot || 0));
      g.imageSmoothingEnabled = true; g.imageSmoothingQuality = 'high';
      g.drawImage(tmp, -M / 2, -M / 2, M, M); g.restore();
    } else g.putImageData(img, 0, 0);
    return g;
  }
  /** 1, 2 or 5 times a power of ten, near a quarter of the window. */
  function niceLen(x) {
    const e = Math.pow(10, Math.floor(Math.log10(x)));
    for (const m of [5, 2, 1]) if (m * e <= x) return m * e;
    return e;
  }
  /** mm -> readable, for scale bars and captions */
  function fmtLen(mm) {
    const a = Math.abs(mm);
    if (!isFinite(mm)) return '--';
    if (a < 1e-10)   return '0';                         // under a milli-angstrom: rounding, not defocus
    if (a >= 1)      return +mm.toPrecision(4) + ' mm';
    if (a >= 1e-3)   return +(mm * 1e3).toPrecision(4) + ' µm';
    if (a >= 1e-6)   return +(mm * 1e6).toPrecision(4) + ' nm';
    return +(mm * 1e7).toPrecision(4) + ' Å';
  }
  function scaleBar(g, N, pxMm, light) {
    const len = niceLen(N * pxMm / 4), w = len / pxMm, s = N / 512;
    g.fillStyle = light ? 'rgba(255,255,255,.75)' : 'rgba(0,0,0,.55)'; g.fillRect(N - w - 26 * s, N - 34 * s, w + 18 * s, 26 * s);
    g.fillStyle = light ? '#1b1f23' : '#fff'; g.fillRect(N - w - 17 * s, N - 16 * s, w, 3 * s);
    g.font = (12 * s) + 'px ui-monospace,Consolas,monospace'; g.textAlign = 'center';
    g.fillText(fmtLen(len), N - w / 2 - 17 * s, N - 20 * s);
  }

  return { wavelengthNm, ABERR, NION_RESIDUALS, cfft2, simulate, simulateDefocused, propagate, propagateField, besselJ, latticeBeams, chooseGrid, buildProbe, specimenPhase,
           fwhmOf, mulberry32, paint, scaleBar, fmtLen, niceLen };
})();

if (typeof window !== 'undefined') window.WaveImages = WaveImages;
