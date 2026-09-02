# STEM Column Ray Trace

Interactive ray diagram of a dedicated STEM column -- condenser, objective, projectors and
EELS coupling, all as **ideal thin lenses** (round lenses only, no aberration terms) -- showing
where the crossovers and diffraction planes land as you change the lens excitations.

**https://seankung.ca/STEM-RayTrace/**

Append `?selftest` to run the built-in assertions against the reference column.

To run it offline, no build step and no npm -- serve the folder with anything static:

```bash
python -m http.server 8123 --directory .
```

## The column

Element positions are **measured off the microscope** and fixed; focal lengths are freely
adjustable (positive only). `z` is measured in mm from the source crossover.

| Element | z (mm) | default f (mm) |
|---|---:|---:|
| Source (point) | 0 | -- |
| C1 | 120 | 120 |
| VOA aperture | 180 | r = 0.030 mm |
| C2 | 240 | 5 |
| C3 | 300 | 55 |
| OL1 | 763 | 22 |
| **Reference plane** (virtual) | 785 | -- |
| OL2 | 807 | 22 |
| PL1 | 910 | 66.6794389630 |
| PL2 | 1030 | off (1/f = 0) |
| PL3 | 1090 | off (1/f = 0) |
| PL4 | 1150 | 41.3326263201 |
| EL (EELS lens) | 1350 | 91.3748845769 |
| **EELS aperture plane** (virtual) | 1550 | -- |
| **EELS focal plane** (virtual) | 1595 | -- |

At the defaults:

* images of the reference plane at **z = 976.68, 1204.28, 1595** -- the EELS focal plane is
  conjugate to it. The illumination crossovers are at **z = 245, 785, 976.68, 1204.28, 1595**,
  so the probe is focused on the reference plane
* diffraction planes at **z = 829, 1181.76, 1550** -- the EELS aperture plane sits on the
  sample's diffraction plane, and 829 is OL2's back focal plane
* PL2 and PL3 are off
* sample to EELS focal plane magnification **1.596x** (= D/L = 45/28.2)
* camera length at the EELS aperture **28.2 mm** (a radian is dimensionless, so r = L*theta
  gives a length; at 15 mrad that puts the direct beam at 0.423 mm radius, which is exactly
  the beam radius printed under the EELS aperture label). This is the value the real column
  runs: at **39 mrad it puts the direct beam at 1.10 mm radius on the aperture**, measured on
  the microscope, and that measurement is what the projector defaults are set from
* convergence semi-angle **alpha = 15 mrad**, set by C2 + C3 against the fixed aperture

The condenser and objective defaults are round, and for a reason: **C1's focal length is the
source distance**, so C1 collimates the source and the VOA sits in a parallel beam where only
its diameter matters; **C2 + C3 sum to their own 60 mm spacing**, so they form an afocal relay;
and **OL1 = OL2 = half the 44 mm objective gap**, so the sample sits at OL1's back focal plane
and OL2's front focal plane. That structure makes the convergence angle exactly

> **alpha = r_VOA * f_C3 / (f_OL1 * f_C2)**, in which f_C1 cancels -- on this column
> **1.364 * (f_C3 / f_C2) mrad**

so the round pair C2 = 5, C3 = 55 gives exactly 15 mrad. The **projectors** get no such luck:
PL1 / PL4 / EL are the exact solve for a 28.2 mm camera length on measured spacings, and those
digits are load-bearing -- rounding them to 4 dp moves the EELS planes by 3e-4 mm.

**Why the objective gap is 44 mm.** The same identity caps the convergence angle at
`alpha_max = r_VOA * 29 / f_OL1 = 0.87 / f_OL1` once C2 hits the 2 mm strength limit. The column
is run at 39 mrad, which needs `f_OL1 <= 22.3 mm`, i.e. a gap of at most about 44.6 mm. At 44 mm
the ceiling is 39.55 mrad and 39 mrad solves with C2 at 2.027 mm, just inside the limit.

A useful consequence of the 45 mm aperture-to-focal gap: **L * M = 45 exactly**, where L is the
camera length at one EELS plane and M the magnification at the other.

(The coupling lens is called **EL**, for EELS lens, leaving CL free to mean condenser lens.)

## How it works

Two engines share one element table, so they cannot drift apart.

* **Ray Optics** (`vendor/rayOptics.js`) traces and renders the real illumination -- the
  ray fan from the point source, clipped by the VOA. Its `IdealLens` is the exact thin
  lens transform, which is what makes the picture trustworthy.
* **A paraxial (ABCD) layer** in `index.html` computes what a ray tracer has no notion of:
  conjugate-plane positions, magnification, camera length, and the auto-focus solves.
  Because ray height at a target plane is exactly affine in a single lens's power, every
  solve is closed form rather than iterative.

### Controls and actions

The lens controls sit **on the figure**, one stack directly above each lens: a **number box** for the
focal length in mm on top, and below it a vertical **slider logarithmic in focal length** --
the number boxes are **staggered into two rows** (OL1 and OL2 are only 44 mm apart, far too close to
put readable boxes side by side), while the sliders stay on one line so the row of thumbs reads as a
power profile across the column -- about
0.86% per step from 10000 mm down to the 2 mm limit, so a 91 mm and a 4 mm lens are equally easy to
dial, with a dedicated **off** (1/f = 0) at the bottom of travel. Up is stronger. Sample defocus is
the horizontal slider under the Sample label, and the alpha slider sits at the top of the Actions
panel. Everything in the
**Display** panel is cosmetic and changes no computed value. The VOA is fixed at 30 um radius.
Under every element label is the **primary beam radius** at that plane, and the readout sits
underneath the figure in four groups.

The **copy** and **save** buttons at the top-right of the figure export it as a PNG at **2x** the
on-screen resolution (about 2100 x 1200 from a 1050 px stage). The controls are HTML layered over
the canvases, so they are never captured in the export. Copying needs a focused window -- if the
browser refuses, the banner says so and save still works.

The action buttons:

| button | solves | for |
|---|---|---|
| Focus probe on sample | C2 + C3 | crossover at the reference plane, at the chosen alpha (alpha = 15 gives C2 5, C3 55) |
| Collimate probe on sample | C2 + C3 | axial slope zero there, at the chosen alpha (alpha = 1 gives C2 3.822, C3 49.831) |
| Match EELS planes | three projector lenses at a time, four if needed | either coupling regime, selected by the toggle above it |

### Reference plane vs specimen

Two distinct planes. The **reference plane** at z = 785, midway between OL1 and OL2, is the
objective's nominal object plane and it never moves: every conjugate, camera length, magnification
and solve is referenced to it, exactly as a real column's projectors are aligned to a fixed height.
The **specimen** sits at 785 + dz and only affects the Bragg cone origin, the beam radius on the
specimen, the probe-defocus readout, and the marker drawn on the figure.

So nothing optical moves when you sweep `dz` -- verified, the EELS solve returns identical
PL1/PL4/EL at dz = 0 and +/- 1 mm. What changes is that the probe is no longer focused *on the
specimen*: the beam there grows as alpha * |dz|, 12 um at dz = 0.8 mm and alpha = 15 mrad.

The `dz` knob spans +/- 1 mm, enough for a confocal depth series.

### Convergence angle

One definition covers both illumination modes:

> **alpha = |marginal ray height at OL2| / f_OL2** -- the slope OL2 imparts to the marginal ray.

Focused, that is the ordinary convergence semi-angle. Parallel, it is the angle OL2 focuses the
illumination to, and the illuminated radius is exactly `r = alpha * f_OL2`, so smaller alpha means a
narrower beam. At the defaults it reads 15 mrad focused; collimating at 1 mrad gives an illuminated
spot of alpha * f_OL2 = 22 um.

The alpha slider is logarithmic over **0.05 - 39.5 mrad** and applies **live**, preserving whichever
mode you last chose with Focus or Collimate. **C2 and C3 set it**, solving both conditions at once
(mode plus angle) in closed form -- the VOA stays fixed at 30 um.

The two modes do **not** have the same reach on this column, and the slider is sized to the focused
one, which is what a dedicated STEM is for:

| mode | reach at f >= 2 mm | what stops it |
|---|---|---|
| focused | **0.047 - 39.5 mrad** | `alpha = 1.364 * f_C3/f_C2` with f_C2 + f_C3 = 60, so the cap gives 1.364*(2/58) at one end and 1.364*(58/2) at the other |
| parallel | **0.0024 - 1.97 mrad** | a 30 um aperture can only be opened out to a 43 um illuminated radius before C2 hits the 2 mm floor |

So **every slider position solves with the probe focused**, and the slider's 0.05 and 39.5 sit just
inside both ends of that range with no dead travel. Above about 1.97 mrad, Collimate refuses and the
banner says why rather than returning something wrong. That asymmetry is the price of a 30 um
aperture: shrinking the aperture is exactly how you reach small convergence angles on a real column,
and it is equally what limits how far the illumination can be spread out. Relaxing the 2 mm cap to
1 mm would move the two ranges to 0.024 - 80.1 and 0.0012 - 3.9 mrad.

The two coupling regimes each necessarily give up the other, so the readout names **which kind** of
plane sits where -- a green tick for the kind the current regime wants, a circle for the other kind,
a parallel sign for collimated illumination -- rather than showing a red cross when you switch mode.
Magnification and camera length are shown at both EELS planes and grey out wherever the plane is not
of the matching kind, since each only means what it is called there.

### Classifying planes

Planes are read off the post-sample transfer matrix [[A,B],[C,D]]: an **image plane** is where
B = 0, a **diffraction plane** is where A = 0. Neither depends on the condenser, so both are correct
in TEM and STEM alike. The **illumination crossovers** in the readout are a different thing --
images of the *source*. With a focused probe they coincide with the sample images; under parallel
illumination they land on the **diffraction** planes instead (829, 1181.76 and 1550 rather than
976.68, 1204.28 and 1595), which is why they must not be used to classify planes.

**Where A and B come from.** Write the ray-transfer matrix from the *reference plane* to whichever
plane you care about as `[[A, B], [C, D]]`. A ray leaving the reference plane with height `y` and
slope `theta` arrives at height `A*y + B*theta`. At a **diffraction plane** `A = 0`, so it arrives at
`B*theta` -- position depends only on angle -- and that `B` *is* the **camera length** `L` in
`r = L*theta`. At an **image plane** `B = 0`, so it arrives at `A*y`, and that `A` *is* the
**magnification**.

Camera length defined this way is the same number as the textbook `L = f_OL2 * M` -- OL2's back focal
plane at z = 829 is conjugate to the aperture with M = -1.2818, so 22 * 1.2818 = 28.2 mm. That is not
a coincidence of the defaults: `A(BFP) = 0` for **any** f_OL2, so writing the BFP -> aperture transfer
as [[a,b],[c,d]] gives `A(aperture) = -b/f_OL2` and `B(aperture) = a*f_OL2 + b*(1 - d1/f_OL2)`. The
diffraction-coupled condition is exactly `b = 0`, which collapses the second to `f_OL2 * a`. So
**`camera length` and `f_OL2 * (magnification of the projector system from OL2's back focal plane)`
are the same quantity**, always, whenever the aperture really is a diffraction plane -- including when
OL2 is detuned away from its default. It is meaningful in STEM as well as TEM: it is what decides
which scattering angles get inside the aperture.

One consequence worth stating on its own, because it is easy to mistake for a geometry statement: with
a focused probe the beam radius printed at the aperture is **exactly `L * alpha`**. It therefore says
nothing about where any element sits -- it is a readout of the camera length you dialled. The shipped
28.2 mm is set from the measured 1.10 mm radius at 39 mrad on the real column.

### The two coupling regimes

They are exact mirrors, and mutually exclusive:

| | diffraction-coupled | image-coupled |
|---|---|---|
| conditions | A(aperture)=0, B(focal)=0 | B(aperture)=0, A(focal)=0 |
| knob | camera length L | magnification M |
| follows | Mag = 45/L at the focal plane | L = 45/M at the focal plane |

Because both conditions live entirely in the post-sample matrix, **matching works identically under
focused and parallel illumination** -- verified, the same PL1 66.68 / PL4 41.33 / EL 91.37 either way.
Verified image-coupled solutions (PL1 / PL4 / EL in mm): M=0.1 -> 60.08 / 6.16 / 98.38 ·
M=1 -> 49.35 / 45.43 / 82.48 · M=5 -> 503.54 / 53.83 / 87.42 · M=10 -> 243.95 / 30.69 / 99.12 ·
M=20 -> 231.09 / 17.15 / 104.43 · M=100 -> 236.48 / 3.87 / 108.92.

### How the solve works

Three conditions need three free lenses, so a third projector lens buys back control of the free
quantity (camera length, or magnification). The solver tries **every triple of projector lenses**
drawn from PL1, PL2, PL3, PL4 and EL, in stages:

0. if the column **already meets** the target, say so and change nothing
1. only lenses that are **already excited** -- so a target that does not need PL2/PL3 never
   disturbs them
2. if nothing works, **switch an off lens on** as a last resort, and say so in the banner
3. still nothing: solve **four at once** (see below)

Every candidate in the winning stage is solved and the **gentlest** one wins, scored by the total
change in lens power. **Hold EL fixed** is on by default and drops EL from the free set, so the EELS lens keeps
whatever value you gave it and the solve uses PL1-PL4 only; combined with the coupling toggle that
gives all four cases.

### Why three lenses, not four

There are exactly three conditions -- the plane condition at the aperture, the plane condition at
the focal plane, and the free quantity (L or M). Three conditions need three unknowns for a locally
unique solution. A fourth lens does **not** over-determine the system; it *under*-determines it into
a one-parameter family.

So "three at a time" does not mean only three lenses are ever used. The solver tries every 3-subset
of the excited projector lenses -- there are five, PL1-PL4 plus EL -- and picks the gentlest, so all
five participate across candidates. What it avoids is varying four *simultaneously*, which has no
unique answer.

That family is still worth having when no triple works, so there is a final stage: scan one lens of
a quad across its range and solve the other three exactly at each step. Measured, it reaches targets
no triple can and leaves EL alone while doing it:

| target | 3 lenses | 4 lenses |
|---|---|---|
| image M = 200 / 800 / 1200 | none | solved |
| diffraction L = 2000 / 8000 / 30000 | none | solved |

The four-lens stage leaves EL alone up to about L = 8000 / M = 2000; past that the all-projector
quads run out and the solver has to bring EL in as well.

It only runs after every triple has failed. Tuned by measurement at 20 outer steps x 1500 inner.
On this column a genuine refusal -- every triple *and* every four-lens family searched -- costs about
**0.13 s with EL held** (the default, one quad) and **0.56 s with EL free** (five quads).

Each triple is tried in **both role orders** -- scanning the downstream lens versus the middle one
gives a different residual with different poles. That swap, not scan density, is what finds the
awkward roots: re-measured on this geometry, a 10x finer root scan recovers nothing the coarse pass
misses (identical reachable sets over M = 0.01..400 and L = 1000..8000) while pushing an exhaustive
three-lens miss from 0.12 s to 0.83 s.

The solve is structured rather than iterative -- alternating three one-dimensional solves does *not*
converge. Both EELS-plane conditions fix the entire ray state at the aperture in closed form,
for diffraction coupling `imgRay = (s*L, -s*L/D)` and `diffRay = (0, -s/L)`, and for image coupling
the mirror `imgRay = (0, s/M)` and `diffRay = (s*M, -s*M/D)`; propagating that back to PL1 and matching the
forward state there makes PL4 an exact affine solve for each trial EL power, leaving one 1-D root
find. Both image parities `s` are searched, which is what makes long camera lengths reachable.
Verified solutions on the default column (PL1 / PL4 / EL in mm):

| L | 1 | 5 | 10 | 20 | 28.2 | 80 | 150 | 400 |
|---|---|---|---|---|---|---|---|---|
| PL1 | 226.70 | 159.26 | 107.16 | 74.18 | **66.68** | 60.22 | 60.00 | 60.24 |
| PL4 | 9.88 | 33.80 | 44.43 | 45.24 | **41.33** | 22.87 | 13.81 | 5.67 |
| EL | 97.40 | 89.84 | 86.31 | 88.19 | **91.37** | 101.67 | 105.42 | 108.31 |
| Mag | 45x | 9x | 4.5x | 2.25x | **1.5957x** | 0.5625x | 0.3x | 0.1125x |

### What bounds the reach

Not the spacing between the projector lenses -- the **lens strength limit**. Measured by sweeping the
cap and probing both regimes with the shipped solver (three-lens stage plus the four-lens fallback):

| strength limit | image-coupled M | diffraction-coupled L | focused alpha |
|---|---|---|---|
| f >= 10 mm | 0.05 - 100 | 0.5 - 1000 mm | 0.27 - 6.7 mrad |
| f >= 5 mm | 0.01 - 800 | 0.05 - 2000 mm | 0.13 - 14.9 mrad |
| **f >= 2 mm** (shipped) | **0.002 - 4000** | **0.01 - 30000 mm** | **0.047 - 39.5 mrad** |
| f >= 1 mm | 0.002 - 8000 | 0.01 - 60000 mm | 0.024 - 80.1 mrad |

Only the innermost part of the EELS range is reachable with three lenses -- at the shipped cap,
**M 0.02-100 and L 0.5-1000 mm**; everything past that needs the four-lens stage, which is why it
exists. 2 mm is stiff but physically defensible, since real magnetic objectives run 1.5-3 mm. On this
column the cap binds on the *condenser* side too, and rather harder: C2 sets both ends of the alpha
slider, and it is the 2 mm floor that makes the 39 mrad working point need a 44 mm objective gap.
Past the ceiling the solver refuses and says so rather than hunting.

### Radial exaggeration

The beam is a fraction of a mm across in a 1595 mm column, so the vertical axis is
magnified by `K`. This is exact rather than a drawing trick: an ideal thin lens is
perfectly linear in (y, dy/dz), so scaling every radial quantity -- ray heights, slopes,
aperture radii -- while leaving z and every focal length untouched reproduces the
identical diagram. Change `K` and no reported plane position moves; `?selftest` asserts it.

### Bragg cones

An optional construction: the illumination re-emitted at +/- theta_B, drawn amber for +1 and violet
for -1. Each cone ray leaves from its own primary ray's position on the sample, so a focused probe
gives cones from a point while parallel illumination gives cones across the whole illuminated width.
The cones merge to a single point at every image plane and spread into separated discs at every
diffraction plane, which is what a CBED pattern shows. They swing considerably wider than the beam,
so switching them on refits the radial scale. They are **on by default at theta_B = 10 mrad**;
since that is below the default alpha = 15 mrad the orders overlap the direct beam, which is the
ordinary overlapping-disc CBED condition. Drop alpha below theta_B and the discs separate instead.
The exact plane positions still come from the paraxial solver, not from the cones.

### Calibrating the camera length against a known d-spacing

The readout's **+/-1 Bragg spacing** is the centre-to-centre distance between the two discs at the
EELS aperture, and it is exactly

> **spacing = 2 * L * theta_B**

because every cone ray is its own primary ray displaced by `L * theta_B`. Three things follow, and
all three are asserted in `?selftest`: the spacing is independent of how wide the direct beam is, it
is the same under focused and parallel illumination, and it does not move with `dz` -- at a
diffraction plane `A = 0`, so `B` is the same measured from the specimen or from the reference plane.

That makes it a **camera-length calibration**. Photograph a known d-spacing, measure the distance
between the +g and -g spots, and invert:

> **L = spacing / (2 * theta_B)**,   equivalently   **L = spacing * d / (2 * lambda)**

**Mind the factor of two at each end.** `theta_B` in this tool is the *deflection* angle of the
scattered beam, `lambda / d` -- twice the Bragg angle. And the direct beam to a single spot is half
the number the readout shows.

Worked example: Si {111}, d = 3.1355 A, at 200 kV (lambda = 2.5079 pm) deflects by
`lambda / d = 8.00 mrad`. Set theta_B = 8.0 and the shipped 28.2 mm camera length predicts
`2 * 28.2 * 0.008 = 0.451 mm` between the +g and -g spots. Measure something else, and the real
camera length is `measured / 0.016` mm.

## Limits

Ideal thin lenses only -- no spherical or chromatic aberration, no space charge, no lens
bore limits, no relativistic or magnetic rotation effects. Positive focal lengths only, and no
lens may be stronger than f = 2 mm.
Element positions are fixed (edit `COLUMN` at the top of the script to change them).

## Licence

The tool is by Sean Kung (<sean.kung@ubc.ca>) and is **MIT** licensed -- see `LICENSE`.

Unlike the other calculators on seankung.ca this one is not a single self-contained file: there is
no build step, but there is one vendored dependency. `vendor/rayOptics.js` is the Ray Optics
Simulation core engine, version 5.4+20260807.c5bf2e2, taken verbatim and unmodified from the
project's `dist-integrations` branch and used under the **Apache License 2.0**. It is committed
rather than fetched from a CDN so the page works offline and the engine version stays pinned.

Apache-2.0 is permissive rather than copyleft: it explicitly permits redistribution, including
publicly and in a work whose own code is MIT, and it does not require this project to adopt the
Apache licence. What it does require is attribution, so travelling with the bundle are:

- `vendor/LICENSE` -- upstream's licence file verbatim. It is the complete Apache-2.0 text, which
  section 4(a) requires be given to recipients, **plus** a `THIRD-PARTY LICENSES` block.
- That block carries the MIT notices for **decimal.js v10.4.3** and **escape-html**, which webpack
  compiled into the bundle. The bundle's banner points at a `rayOptics.js.LICENSE.txt` sidecar that
  upstream does not publish on this branch; those notices live in `vendor/LICENSE` instead.
- `NOTICE` -- the summary of all of the above.

Nothing here is legal advice, but the obligations are the plain reading of the licence text.
