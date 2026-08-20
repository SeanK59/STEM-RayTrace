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

Element positions are fixed; focal lengths are freely adjustable (positive only).
`z` is measured in mm from the source crossover.

| Element | z (mm) | default f (mm) |
|---|---:|---:|
| Source (point) | 0 | -- |
| C1 | 100 | 100 |
| VOA aperture | 200 | r = 0.4 mm |
| C2 | 300 | 100 |
| C3 | 450 | 50 |
| OL1 | 600 | 40 |
| **Reference plane** (virtual) | 640 | -- |
| OL2 | 680 | 40 |
| PL1 | 800 | 120 |
| PL2 | 900 | off (1/f = 0) |
| PL3 | 1000 | off (1/f = 0) |
| PL4 | 1100 | 150 |
| EL (EELS lens) | 1250 | 375 |
| **EELS aperture plane** (virtual) | 1300 | -- |
| **EELS focal plane** (virtual) | 1500 | -- |

At the defaults:

* images of the reference plane at **z = 920, 1500** -- the EELS focal plane is conjugate
  to it. The illumination crossovers are at **z = 400, 640, 920, 1500**, so the probe is
  focused on the reference plane
* diffraction planes at **z = 720, 1300** -- the EELS aperture plane sits on the
  sample's diffraction plane
* PL2 and PL3 are off
* sample to EELS focal plane magnification **5.0x**
* camera length at the EELS aperture **40 mm** (a radian is dimensionless, so r = L*theta
  gives a length; at 5 mrad that puts the direct beam at 0.2 mm radius, which is exactly
  the beam radius printed under the EELS aperture label)
* convergence semi-angle **alpha = 5 mrad**, set by C2 + C3 against the fixed aperture

The numbers are round because the structure is physically natural: C1 collimates the
source (so the VOA sits at C1's back focal plane), C2 and C3 relay and re-collimate,
the sample sits at OL1's back focal plane and OL2's front focal plane, PL4 collimates
the axial ray, and EL forms the sample image at its own back focal plane.

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
focal length in mm on top, and below it a vertical **slider logarithmic in focal length** -- about
0.86% per step from 10000 mm down to the 2 mm limit, so a 375 mm and a 4 mm lens are equally easy to
dial, with a dedicated **off** (1/f = 0) at the bottom of travel. Up is stronger. Sample defocus is
the horizontal slider under the Sample label, and the alpha slider sits at the top of the Actions
panel. Everything in the
**Display** panel is cosmetic and changes no computed value. The VOA is fixed at 0.40 mm radius.
Under every element label is the **primary beam radius** at that plane, and the readout sits
underneath the figure in four groups.

The **copy** and **save** buttons at the top-right of the figure export it as a PNG at **2x** the
on-screen resolution (about 2100 x 1200 from a 1050 px stage). The controls are HTML layered over
the canvases, so they are never captured in the export. Copying needs a focused window -- if the
browser refuses, the banner says so and save still works.

The action buttons:

| button | solves | for |
|---|---|---|
| Focus probe on sample | C2 + C3 | crossover at the reference plane, at the chosen alpha (alpha = 5 gives C2 100, C3 50) |
| Collimate probe on sample | C2 + C3 | axial slope zero there, at the chosen alpha (alpha = 20/11 gives C2 100, C3 34.375) |
| Match EELS planes | three projector lenses at a time, four if needed | either coupling regime, selected by the toggle above it |

### Reference plane vs specimen

Two distinct planes. The **reference plane** at z = 640, midway between OL1 and OL2, is the
objective's nominal object plane and it never moves: every conjugate, camera length, magnification
and solve is referenced to it, exactly as a real column's projectors are aligned to a fixed height.
The **specimen** sits at 640 + dz and only affects the Bragg cone origin, the beam radius on the
specimen, the probe-defocus readout, and the marker drawn on the figure.

So nothing optical moves when you sweep `dz` -- verified, the EELS solve returns identical
PL1/PL4/EL at dz = 0 and +/- 1 mm. What changes is that the probe is no longer focused *on the
specimen*: the beam there grows as alpha * |dz|, 4 um at dz = 0.8 mm and alpha = 5 mrad.

The `dz` knob spans +/- 1 mm, enough for a confocal depth series.

### Convergence angle

One definition covers both illumination modes:

> **alpha = |marginal ray height at OL2| / f_OL2** -- the slope OL2 imparts to the marginal ray.

Focused, that is the ordinary convergence semi-angle. Parallel, it is the angle OL2 focuses the
illumination to, and the illuminated radius is exactly `r = alpha * f_OL2`, so smaller alpha means a
narrower beam. At the defaults it reads 5 mrad focused; collimating without touching alpha keeps
it at 5 mrad and widens the illuminated spot to alpha * f_OL2 = 0.2 mm.

The alpha slider is logarithmic over **0.55 - 50 mrad** and applies **live**, preserving whichever
mode you last chose with Focus or Collimate. **C2 and C3 set it**, solving both conditions at once
(mode plus angle) in closed form -- the VOA stays fixed at 0.4 mm. With that aperture the focused
floor is 0.137 mrad, comfortably below the slider's 0.55 mrad minimum, so every slider position
solves. (At the 1.6 mm aperture this tool shipped with originally the floor was 0.549 mrad, right at
the bottom of the slider -- shrinking the aperture is exactly how you reach smaller angles on a real
column, which is why alpha = 5 mrad is set that way here rather than by straining the condensers.)

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
illumination they land on the **diffraction** planes instead (720 and 1300 rather than 920 and
1500), which is why they must not be used to classify planes.

**Where A and B come from.** Write the ray-transfer matrix from the *reference plane* to whichever
plane you care about as `[[A, B], [C, D]]`. A ray leaving the reference plane with height `y` and
slope `theta` arrives at height `A*y + B*theta`. At a **diffraction plane** `A = 0`, so it arrives at
`B*theta` -- position depends only on angle -- and that `B` *is* the **camera length** `L` in
`r = L*theta`. At an **image plane** `B = 0`, so it arrives at `A*y`, and that `A` *is* the
**magnification**.

Camera length defined this way is the same number as the textbook `L = f_OL2 * M` -- OL2's back focal
plane at z = 720 is conjugate to the aperture with M = -1, so 40 * 1 = 40 mm. It is meaningful in
STEM as well as TEM: it is what decides which scattering angles get inside the aperture.

### The two coupling regimes

They are exact mirrors, and mutually exclusive:

| | diffraction-coupled | image-coupled |
|---|---|---|
| conditions | A(aperture)=0, B(focal)=0 | B(aperture)=0, A(focal)=0 |
| knob | camera length L | magnification M |
| follows | Mag = 200/L at the focal plane | L = 200/M at the focal plane |

Because both conditions live entirely in the post-sample matrix, **matching works identically under
focused and parallel illumination** -- verified, the same PL1 120 / PL4 150 / EL 375 either way.
Verified image-coupled solutions (PL1 / PL4 / EL in mm): M=1 -> 26.09 / 90 / 12.10 ·
M=10 -> 243.98 / 47.87 / 69.34 · M=20 -> 265.04 / 31.14 / 79.39.

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
| image M = 800 / 1200 / 5000 | none | solved |
| diffraction L = 16000 / 30000 | none | solved |

It only runs after every triple has failed. Tuned by measurement at 20 outer steps x 1500 inner:
that keeps every target a richer 25 x 3000 setting reaches, while a genuine refusal costs 0.52 s
with EL held (the default, one quad) and 2.3 s with EL free (five quads).

Each triple is tried in **both role orders** -- scanning the downstream lens versus the middle one
gives a different residual with different poles. That swap, not scan density, is what finds the
awkward roots: a 10x finer scan was measured to recover nothing the coarse pass misses while
pushing an exhaustive miss from 0.19 s to 1.9 s.

The solve is structured rather than iterative -- alternating three one-dimensional solves does *not*
converge. Both EELS-plane conditions fix the entire ray state at the aperture in closed form,
for diffraction coupling `imgRay = (s*L, -s*L/D)` and `diffRay = (0, -s/L)`, and for image coupling
the mirror `imgRay = (0, s/M)` and `diffRay = (s*M, -s*M/D)`; propagating that back to PL1 and matching the
forward state there makes PL4 an exact affine solve for each trial EL power, leaving one 1-D root
find. Both image parities `s` are searched, which is what makes long camera lengths reachable.
Verified solutions on the default column (PL1 / PL4 / EL in mm):

| L | 10 | 20 | 40 | 80 | 150 | 400 |
|---|---|---|---|---|---|---|
| PL1 | 400 | 240 | **120** | 81.8 | 70.8 | 61.5 |
| PL4 | 200 | 150 | **150** | 200 | 729.7 | 30.9 |
| EL | 62.5 | 125 | **375** | 281.3 | 164.9 | 78.5 |
| Mag | 20x | 10x | **5x** | 2.5x | 1.33x | 0.5x |

### What bounds the reach

Not the spacing between the projector lenses -- the **lens strength limit**. Reach scales as roughly
**M_max ~ 800 / f_min**, measured by sweeping the cap:

| strength limit | image-coupled M | diffraction-coupled L |
|---|---|---|
| f >= 10 mm | 0.05 - 80 | 0.5 - 1200 mm |
| f >= 5 mm | 0.02 - 120 | 0.5 - 2000 mm |
| **f >= 2 mm** (shipped) | **0.02 - 400** | **0.5 - 8000 mm** |
| f >= 1 mm | 0.02 - 800 | 0.5 - 8000+ mm |

At any cap the working lens sits hard against the floor while the others barely move -- at f >= 10 mm
it was PL4 at 10.12 mm for M = 80 and 13.04 mm for L = 1200, with PL1 near 290 mm and EL near 90 mm
either way. 2 mm is stiff but physically defensible, since real magnetic objectives run 1.5-3 mm.
Past the ceiling the solver refuses and says so rather than hunting.

### Radial exaggeration

The beam is a fraction of a mm across in a 1500 mm column, so the vertical axis is
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
since that is below alpha = 20 mrad the orders overlap the direct beam, which is the ordinary
overlapping-disc CBED condition. The exact plane positions still come from the paraxial
solver, not from the cones.

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
