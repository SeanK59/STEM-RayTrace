# STEM Column Ray Trace

Interactive ray diagram of a dedicated STEM column -- condenser, objective, projectors and
EELS coupling, all as **ideal thin lenses** (round lenses only, no aberration terms) -- showing
where the crossovers and diffraction planes land as you change the lens excitations.

**https://seankung.ca/STEM-RayTrace/**

The page itself carries only a short "how to read it"; the usage notes and the optics behind the controls
live in `help.html`, linked from the header and footer.

Append `?selftest` to run the built-in assertions against the reference column.

To run it offline, no build step and no npm -- serve the folder with anything static:

```bash
python -m http.server 8123 --directory .
```

## The column

Element positions are **measured off the microscope** and fixed; focal lengths are freely
adjustable (positive only). `z` is measured in mm from the source crossover. The table lives in
`column.js`, the one file both pages of this repo load.

| Element | z (mm) | default f (mm) |
|---|---:|---:|
| Source (point) | 0 | -- |
| C1 | 120 | off (1/f = 0) |
| VOA aperture | 180 | 30 um diameter (r = 0.015 mm) |
| C2 | 240 | 21.8181818182 (= 240/11) |
| C3 | 300 | 36 |
| OL1 (objective pre-field) | 783 | 2 |
| **Reference plane** (virtual) | 785 | -- |
| OL2 (objective post-field) | 787 | 2 |
| PL1 | 910 | 56.0828398219 |
| PL2 | 1030 | 13.5835644681 |
| PL3 | 1090 | off (1/f = 0) |
| PL4 | 1150 | 14.0309991241 |
| EL (EELS lens) | 1350 | 105 |
| **EELS aperture plane** (virtual) | 1550 | -- |
| **EELS focal plane** (virtual) | 1595 | -- |

No round lens may be stronger than **f = 4 mm**; the objective's limit is **1.5 mm**.

At the defaults:

* images of the reference plane at **z = 966.08, 1047.25, 1166.25, 1595** -- the EELS focal plane is
  conjugate to it. The illumination crossovers are at **z = 264, 785**, then the same image planes,
  so the probe is focused on the reference plane
* diffraction planes at **z = 789, 1014.53, 1141.58, 1550** -- the EELS aperture plane sits on the
  sample's diffraction plane, and 789 is OL2's back focal plane
* C1 and PL3 are off
* sample to EELS focal plane magnification **1.596x** (= D/L = 45/28.2)
* camera length at the EELS aperture **28.2 mm** (a radian is dimensionless, so r = L*theta gives a
  length; at 15 mrad that puts the direct beam at 0.423 mm radius, which is exactly the beam radius
  printed under the EELS aperture label). This is the value the real column runs: at **39 mrad it
  puts the direct beam at 1.10 mm radius on the aperture**, measured on the microscope, and that
  measurement is what the projector defaults are set from
* convergence semi-angle **alpha = 15 mrad**, set by C2 + C3 against the fixed aperture with C1 off:
  C2 makes a crossover at 264, which is C3's front focal point, so C3 sends a parallel beam to OL1

### The objective

The UltraSTEM 200's Riecke-Ruska condenser-objective: a **4 mm pole-piece gap** and **f = 2 mm** for
both its pre-field and its post-field (Dellby et al. 2011). It is modelled as two thin lenses on the
pole faces, OL1 and OL2, 2 mm either side of the specimen, and driven by **one control**: its halves are
always equally excited, and 4 mm is too close for two control stacks anyway.

With f = 2 mm on the pole faces, the specimen sits **exactly on OL1's back focal plane** and on OL2's
front focal plane. A parallel beam into OL1 is focused on the specimen, so the convergence angle is the
beam radius at OL1 over f_OL1. OL2 sends the specimen's image to infinity -- the first real image is
PL1's -- and puts the diffraction pattern on its back focal plane at 789. Parallel illumination needs
the condenser crossover on OL1's front focal plane, at 781.

A condenser beam that is not quite parallel moves the probe crossover off the specimen by Newton's
f^2/x, for a condenser crossover a distance x from OL1's front focal plane -- about 8 um for one near
C3, on one side or the other depending on which side of that plane it sits. That is the scale of the
"sample-level beam crossover a few um above the center of the objective lens" that Dellby et al.
describe as normal operation.

Krivanek et al. 2008 quote f = 1.5 mm, a 1.5 mm bore and Cc about 1 mm for the earlier 100 kV
objective; Dellby et al. describe the lens redesigned for 200 kV. The **VOA** is 30 um in diameter,
the value Krivanek et al. give as typical.

The **projectors** are the exact solve for a 28.2 mm camera length with EL at a round 105 mm. EL is
the free parameter of a one-parameter family of four-lens solutions, and 105 mm is its gentlest
member -- strongest lens 13.58 mm. The digits are load-bearing: rounding them to 4 dp moves the EELS
planes by about 1e-4 mm.

A useful consequence of the 45 mm aperture-to-focal gap: **L * M = 45 exactly**, where L is the
camera length at one EELS plane and M the magnification at the other.

(The coupling lens is called **EL**, for EELS lens, leaving CL free to mean condenser lens. On the
real column that job is done by a quadrupole-octupole coupling module; see Limits.)

References:

- N. Dellby et al., "Dedicated STEM for 200 to 40 keV operation", *Eur. Phys. J. Appl. Phys.* **54**,
  33505 (2011).
- O. L. Krivanek et al., "An electron microscope for the aberration-corrected era",
  *Ultramicroscopy* **108**, 179-195 (2008).

## How it works

Two engines share one element table, so they cannot drift apart.

* **Ray Optics** (`vendor/rayOptics.js`) traces and renders the real illumination -- the
  ray fan from the point source, clipped by the VOA. Its `IdealLens` is the exact thin
  lens transform, which is what makes the picture trustworthy.
* **A paraxial (ABCD) layer** in `index.html` computes what a ray tracer has no notion of:
  conjugate-plane positions, magnification, camera length, and the auto-focus solves.
  Because ray height at a target plane is exactly affine in a single lens's power, every
  solve is closed form rather than iterative.

`column.js` holds the element table, the strength limits, and a `transfer()` that returns the
paraxial matrix between any two planes. The self-test asserts it agrees with the page's own
tracer element for element.

### Controls and actions

The lens controls sit **on the figure**, one stack directly above each lens: a **number box** for the
focal length in mm on top, and below it a vertical **slider logarithmic in focal length**. The number
boxes are **staggered into two rows**, since C2/C3 and PL2/PL3/PL4 are only 60 mm apart, too close for
readable boxes side by side. The sliders stay on one line, so the row of thumbs reads as a power
profile across the column. They run at about 0.8% per step from 10000 mm down to the strength limit,
so a 105 mm and a 4 mm lens are equally easy to dial, with a dedicated **off** (1/f = 0) at the bottom of
travel. Up is stronger. The objective's single stack, labelled **OL**, sits at the specimen between its
two halves and runs down to 1.5 mm (it normally sits at 2 mm). Sample defocus is the horizontal slider under the Sample label, and
the alpha slider sits at the top of the Actions panel. The **Display** panel is cosmetic with one
exception: theta_B sets the +/-1 Bragg spacing reported at the aperture. Nothing in it moves a plane. The
VOA is fixed at 30 um diameter. Under every element label is the **primary beam radius** at that plane, and
the readout sits underneath the figure in four groups.

The **checkbox at the foot of a lens stack pins that lens**: a solve may not move a pinned lens.
**EL is pinned on a fresh load** -- holding the EELS lens while PL1-PL4 do the work is the common
case -- and nothing else is. Only the eight lenses a solve can reach carry one: C1, C2 and C3, owned
by the illumination solve, and PL1-PL4 plus EL, owned by the EELS solve. The objective is set by hand
only, so a pin there would mean nothing and it has none. A pinned lens is drawn with an accent border
and a dimmed slider, so the pinned set reads at a glance.

**Horizontal zoom** sits beside copy/save: `-`, a log slider, `+`, a readout, **reset** (back to the
whole column) and **fit K**. It stretches the z axis only, from 1x to 2000x -- at the top of the range the
window is 0.86 mm wide, so the 4 mm between OL1 and OL2 can be filled with beam -- and the radial
exaggeration `K` is not touched; no plane, camera length or other reported number moves. A slim scroll bar
under the figure, dragging the figure, sideways scroll and shift + wheel pan along the column, and Ctrl/Cmd
+ wheel (or a trackpad pinch) zooms about the cursor. Because the beam is only ~30 um across at the
objective, **fit K** refits `K` to the beam inside the current window (by click only, never while
panning). The canvases stay viewport-sized, so there is no huge bitmap at high zoom.

Two details make this work with the ray engine, which has one scale for both axes and draws its line
widths in scene units: the radial exaggeration rides in the engine's y units as `K / zoom`, so a y unit
stays `K * sy` px tall however far z is stretched, and `lengthScale` is set to `1 / zoom` so rays stay
~1 px thick instead of growing with the zoom.

The **copy** and **save** buttons at the top-right of the figure export it as a PNG at **2x** the
on-screen resolution (about 2100 x 1200 from a 1050 px stage), showing whatever window is on screen. The controls are HTML layered over
the canvases, so they are never captured in the export. Copying needs a focused window -- if the
browser refuses, the banner says so and save still works.

The action buttons:

| button | solves | for |
|---|---|---|
| Focus probe on sample | C2 + C3, and C1 when they cannot | crossover at the reference plane, at the chosen alpha (alpha = 15 gives C2 240/11, C3 36, C1 off) |
| Collimate probe on sample | C2 + C3, and C1 when they cannot | axial slope zero there, at the chosen alpha (alpha = 0.3 gives C2 7.088, C3 47.493, C1 off) |
| **Solve** | three projector lenses at a time, four if needed -- or **two** with no target | matches both EELS planes, in the coupling regime selected by the toggle above it |

### Reference plane vs specimen

Two distinct planes. The **reference plane** at z = 785, midway between OL1 and OL2, is the
objective's nominal object plane and it never moves: every conjugate, camera length, magnification
and solve is referenced to it, exactly as a real column's projectors are aligned to a fixed height.
The **specimen** sits at 785 + dz and only affects the Bragg cone origin, the beam radius on the
specimen, the probe-defocus readout, and the marker drawn on the figure.

So nothing optical moves when you sweep `dz` -- verified, the EELS solve returns identical
projector settings at dz = 0 and +/- 1 mm. What changes is that the probe is no longer focused *on
the specimen*: the beam there grows as alpha * |dz|, 15 um at dz = 1 mm and alpha = 15 mrad.

The `dz` knob spans +/- 1 mm -- inside the 4 mm gap, and enough for a confocal depth series.

### Convergence angle

> **alpha = the marginal ray's slope inside the objective gap.**

No lens sits between OL1 and OL2, so that slope is the same at the reference plane and on the
specimen, and with the probe focused it is the ordinary convergence semi-angle -- the beam radius at
OL1 over f_OL1, since the beam arrives there parallel. Under parallel illumination that slope is zero,
so there alpha means the angle OL2 focuses the beam to instead, and the illuminated radius is exactly
`r = alpha * f_OL2`; collimating at 0.3 mrad gives a 0.6 um spot.

The alpha slider is logarithmic over **0.05 - 50 mrad** and applies **live**, preserving whichever
mode you last chose with Focus or Collimate. **C2 and C3 set it**, solving both conditions at once
(mode plus angle) in closed form -- the VOA stays fixed. When they cannot, the solve **brings
C1 in**: it scans C1 across its range, solves C2 + C3 exactly at every step, and keeps the setting
whose strongest condenser lens is weakest. The banner says when it does. A C1 that is already on stays
where it is whenever C2 + C3 can follow, so C1 only moves when it has to.

| mode | C2 + C3 alone (C1 off) | with C1 free |
|---|---|---|
| focused | **0.714 - 137.5 mrad** | from **0.043 mrad** upward, no ceiling found below 1000 mrad |
| parallel | **0.0030 - 0.57 mrad** | from **0.00028 mrad** upward, to at least 100 mrad |

So **every slider position solves in both modes**, and the column's **39 mrad working point needs no
C1** -- C2 11.65, C3 47.76. The edges of the C2 + C3 window are C2 or C3 reaching the 4 mm floor;
outside it C1 comes on. Pin C1 and that second stage is off: outside the window the solve then refuses
and says why.

Every candidate is also **traced forward and checked** before it is accepted. The closed form is exact
in principle, but a C1 crossover a few mm ahead of C2 chains several large magnifications together,
and a matching slip of 1e-6 mm there once turned into a visibly wrong probe.

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
illumination they land on the **diffraction** planes instead (789, 1014.53, 1141.58 and 1550
rather than 966.08, 1047.25, 1166.25 and 1595), which is why they must not be used to classify planes.

**Where A and B come from.** Write the ray-transfer matrix from the *reference plane* to whichever
plane you care about as `[[A, B], [C, D]]`. A ray leaving the reference plane with height `y` and
slope `theta` arrives at height `A*y + B*theta`. At a **diffraction plane** `A = 0`, so it arrives at
`B*theta` -- position depends only on angle -- and that `B` *is* the **camera length** `L` in
`r = L*theta`. At an **image plane** `B = 0`, so it arrives at `A*y`, and that `A` *is* the
**magnification**.

Camera length defined this way is the same number as the textbook `L = f_OL2 * M` -- OL2's back focal
plane at z = 789 is conjugate to the aperture with M = -14.1, so 2 * 14.1 = 28.2 mm. That is not a
coincidence of the defaults: `A(BFP) = 0` for **any** f_OL2 and any specimen-to-OL2 distance, so
writing the BFP -> aperture transfer as [[a,b],[c,d]] gives `A(aperture) = -b/f_OL2`, and the
diffraction-coupled condition is exactly `b = 0`, which collapses `B(aperture)` to `f_OL2 * a`. So
**`camera length` and `f_OL2 * (magnification of the projector system from OL2's back focal plane)`
are the same quantity**, always, whenever the aperture really is a diffraction plane. It is meaningful
in STEM as well as TEM: it is what decides which scattering angles get inside the aperture. It is also
why long camera lengths are the hard end of the EELS reach on this column: with f_OL2 = 2 mm, every
millimetre of camera length has to be bought as projector magnification.

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
focused and parallel illumination** -- verified, the same PL1 56.08 / PL2 13.58 / PL4 14.03 / EL 105
either way. Verified image-coupled solutions from the shipped column with EL free (mm):

| M | lenses used | PL1 | PL2 | PL4 | EL |
|---:|---|---:|---:|---:|---:|
| 1 | PL2 + PL4 + EL | 56.08 | 13.47 | 11.44 | 96.68 |
| 5 | PL1 + PL4 + EL | 55.22 | 13.58 | 34.51 | 85.07 |
| 10 | PL1 + PL4 + EL | 54.31 | 13.58 | 127.56 | 157.06 |
| 20 | PL1 + PL4 + EL | 43.30 | 13.58 | 59.12 | 47.38 |
| 50 | PL1 + PL2 + EL | 24.98 | 50.82 | 14.03 | 110.33 |
| 100 | PL1 + PL2 + EL | 51.89 | 42.26 | 14.03 | 105.40 |
| 200 | PL1 + PL2 + EL | 73.69 | 32.73 | 14.03 | 104.96 |

### How the solve works

Three conditions need three free lenses, so a third projector lens buys back control of the free
quantity (camera length, or magnification). The solver tries **every triple of projector lenses**
drawn from PL1, PL2, PL3, PL4 and EL, in stages:

0. if the column **already meets** the target, say so and change nothing
1. only lenses that are **already excited** -- so a target that does not need PL3 never disturbs it
2. if nothing works, **switch an off lens on** as a last resort, and say so in the banner
3. still nothing: solve **four at once** (see below)

Every candidate in the winning stage is solved and the **gentlest** one wins, scored by the total
change in lens power.

**Pinning** removes a lens from the free set for the whole search, including the four-lens stage, so
the solve works around whatever you have fixed. Three conditions need three unknowns, so the EELS
solve needs **at least three free projector lenses**; pin more than two of the five and it refuses
up front and names what is pinned rather than grinding through a search that cannot succeed. The
**Hold EL fixed** checkbox in the Actions panel is the same pin as EL's own box on the figure --
one flag, two views -- kept because holding the EELS lens is the common case; combined with the
coupling toggle it gives all four cases.

The illumination solve has no slack in C2 and C3: it sets the mode *and* the convergence angle with
exactly those two, so pinning either one stops it outright with a message saying so. Pinning C1
only stops it from bringing C1 in.

### Solving without a target

The checkbox beside the camera length decides whether that value is an input at all. It is **on by
default**, which is the behaviour above: L (or M) is a third condition, and three conditions need
three free projector lenses.

Clear it and the target stops being a condition. The solve then matches the two planes on their own,
which is **two conditions and therefore two lenses** -- every pair, in both role orders, gentlest
wins -- and whatever camera length falls out is **written back to the slider**, so the control turns
into a readout. Turning the checkbox back on then aims at what was just found.

Whether a given pair has a solution depends on where the lenses you did *not* free are sitting. From
the shipped column, with EL held, the other two projectors held at their defaults, and the first lens
of the pair knocked to 40 mm to start:

| free pair | L (mm) | the two focal lengths |
|---|---:|---|
| PL1 + PL4 | **28.20** | 56.083 / 14.031 -- the shipped column itself |
| PL2 + PL4 | **28.20** | 13.584 / 14.031 -- likewise |
| PL2 + PL3 | 5.26 | 24.990 / 34.096 |
| PL1 + PL2 | 3.66 | 140.838 / 100.165 |
| PL1 + PL3 | 0.68 | 100.431 / 24.536 |
| PL3 + PL4 | none | -- |

With the other two projectors switched off instead, a different three solve -- PL1 + PL4
(L = 8.39 mm), PL2 + PL4 (2.35 mm) and PL3 + PL4 (1.25 mm) -- so two is the degree-of-freedom count,
not a guarantee. The four-lens stage does
not run without a target: it exists to extend *reach* toward a value you named, and there is no value
to reach.

### Why three lenses, not four

There are exactly three conditions -- the plane condition at the aperture, the plane condition at
the focal plane, and the free quantity (L or M). Three conditions need three unknowns for a locally
unique solution. A fourth lens does **not** over-determine the system; it *under*-determines it into
a one-parameter family. The shipped projector setting is itself a member of such a family.

So "three at a time" does not mean only three lenses are ever used. The solver tries every 3-subset
of the excited projector lenses -- there are five, PL1-PL4 plus EL -- and picks the gentlest, so all
five participate across candidates. What it avoids is varying four *simultaneously*, which has no
unique answer.

That family is still worth having when no triple works, so there is a final stage: scan one lens of
a quad across its range and solve the other three exactly at each step. Measured from the shipped
column with EL free, it reaches targets no triple can:

| target | 3 lenses | 4 lenses |
|---|---|---|
| image M = 5000 / 10000 / 20000 | none | solved, EL untouched |
| diffraction L = 200 / 400 / 600 | none | solved, EL untouched |
| diffraction L = 0.002 - 0.007, image M = 0.07 - 0.2 | none | solved, moving EL |

It only runs after every triple has failed. Tuned by measurement at 20 outer steps x 1500 inner: a
richer 25 x 3000 reaches exactly the same targets at about twice the cost. On this column a genuine
refusal -- every triple *and* every four-lens family searched -- costs about **0.3 - 0.4 s with EL
held** (the default, one quad) and **1.6 - 2.4 s with EL free** (five quads).

Each triple is tried in **both role orders** -- scanning the downstream lens versus the middle one
gives a different residual with different poles. That swap, not scan density, is what finds the
awkward roots: re-measured on this geometry, a 10x finer root scan recovers nothing the coarse pass
misses (L = 800, M = 0.05 and M = 40000 all still miss) while pushing an exhaustive miss from about
1.8 s to 4 - 5 s.

The solve is structured rather than iterative -- alternating three one-dimensional solves does *not*
converge. Both EELS-plane conditions fix the entire ray state at the aperture in closed form,
for diffraction coupling `imgRay = (s*L, -s*L/D)` and `diffRay = (0, -s/L)`, and for image coupling
the mirror `imgRay = (0, s/M)` and `diffRay = (s*M, -s*M/D)`; propagating that back to the most
upstream free lens and matching the forward state there makes the middle lens an exact affine solve
for each trial power of the downstream one, leaving one 1-D root find. Both image parities `s` are
searched, which is what makes long camera lengths reachable. Verified solutions from the shipped
column with EL free (mm):

| L | 1 | 5 | 10 | 20 | 28.2 | 50 | 80 |
|---|---|---|---|---|---|---|---|
| lenses used | PL1+PL2+EL | PL1+PL4+EL | PL1+PL4+EL | PL1+PL4+EL | **shipped** | PL1+PL2+EL | PL1+PL2+EL |
| PL1 | 35.78 | 56.90 | 56.21 | 56.09 | **56.08** | 57.81 | 58.70 |
| PL2 | 44.97 | 13.58 | 13.58 | 13.58 | **13.58** | 8.58 | 5.69 |
| PL4 | 14.03 | 39.87 | 28.67 | 18.24 | **14.03** | 14.03 | 14.03 |
| EL | 98.67 | 86.97 | 96.88 | 103.04 | **105** | 105.06 | 105.08 |
| Mag | 45x | 9x | 4.5x | 2.25x | **1.5957x** | 0.9x | 0.5625x |

### What bounds the reach

Not the spacing between the projector lenses -- the **lens strength limit**, and behind it the strong
objective. Measured with the shipped solver from the shipped column, every stage included:

| | three lenses | all stages |
|---|---|---|
| diffraction-coupled L, EL free | 0.02 - 100 mm | **0.002 - 700 mm** |
| diffraction-coupled L, EL held (default) | 2 - 100 mm | **0.05 - 600 mm** |
| image-coupled M, EL free | 0.5 - 3000 | **0.06 - 20000** |
| image-coupled M, EL held (default) | 20 - 3000 | **2 - 10000** |

Past those ends the solver refuses and says so rather than hunting. The ends are sampled, so read them
as "about".

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

## Wave-optics view

`wave.html` is the same column -- the same `column.js`, the same thin lenses, the same lens settings --
with electron-wave mechanics added and nothing else. Open it from the ray page's **wave-optics view**
link and it starts from that page's lenses (they travel in the URL hash); leave both open and it
follows every change live over a `BroadcastChannel`. The ray page stays the control panel. The wave
page adds only what a ray cannot carry: the accelerating voltage (30 / 80 / 200 kV), the probe's
residual aberrations, and a specimen.

It draws two planes:

* **the specimen** -- the probe intensity, or the specimen's phase
* **any plane downstream** -- the EELS aperture (default), the EELS focal plane, or any z

**The probe** is the converging wave the ray model describes, a sharp cone of half-angle alpha, with
the aberration phase in Krivanek notation up to third order:
`chi = (2 pi / lambda) * sum theta^(n+1)/(n+1) * [C_nm,a cos(m phi) + C_nm,b sin(m phi)]` for
C10, C12, C21, C23, C30, C32, C34, in nm. The **Nion residuals** button loads 2 nm at first order,
20 nm at second and 200 nm at third. The ray page's probe defocus -- how far the specimen sits past the
crossover -- is added to C10, so moving the specimen on the ray page defocuses the probe here. Under
parallel illumination the incident wave is a plane wave.

**Propagation** uses the column's own transfer matrix from the specimen, through the Collins
integral, which is the wave-optical form of the ABCD matrix the ray page draws with. It is evaluated
two ways, each exact where it samples well, split at `|B/A| = z* = N dx^2 / lambda`:

* near a diffraction plane: a pre-chirp and one FFT (at A = 0, a pure scaled Fourier transform)
* near an image plane: an angular-spectrum propagation over B/A, then a magnification A

At the EELS aperture in diffraction coupling that gives CBED discs of radius L*alpha at spacing
L*lambda/d, exactly where the ray page puts the beam edge and its Bragg discs (the optional dashed
overlay). The EELS limiting aperture is ignored: every electron reaches every plane.

**The window** is chosen per render so the pixel resolves every angle present -- the probe cone plus
two Bragg orders, times an oversample factor -- and is then N pixels of that, nudged to a whole number
of lattice periods. A probe defocused past what the window holds is flagged and not drawn; at 20 mrad
and 512 pixels that is about 250 nm of defocus, far less than the ray page's +/- 1 mm `dz`.
Use C10 for fine defocus, or raise the grid. A 512 x 512 render takes about 70 - 90 ms, and
1024 x 1024 about 0.4 s.

**Specimens** are thin phase objects:

* a square lattice, `t = exp(i sigma [cos(2 pi x/d) + cos(2 pi y/d)])`, whose orders sit at
  +/- lambda/d along x and y (Si {111}, d = 3.1355 A, is the default)
* amorphous carbon, as a Gaussian-correlated random phase screen
* vacuum

There is no multislice, absorption, inelastic scattering, partial coherence or chromatic spread.

`wave.html?selftest` checks the machinery against closed-form results:

* the FFT (Parseval, and a round trip to 4e-16)
* the aberration-free probe's FWHM against the Airy value 0.5145 lambda/alpha (within 0.1%)
* a Gaussian beam propagated to five planes, spanning both evaluations, against its analytic width
  (agreement to better than 1e-6)
* the CBED disc radius against L*alpha, and the +/-1 spacing against 2 L lambda/d (each within one
  pixel)
* the 1.10 mm direct-beam radius at 39 mrad
* the out-of-window flag
* the hand-off decoding

The ray page's own self-test asserts that `column.js`'s `illumination()` and `transfer()` agree with
its tracer in every regime, so the two pages share one geometry by construction.

## Limits

Ideal thin lenses only -- no spherical or chromatic aberration, no space charge, no lens
bore limits, no relativistic or magnetic rotation effects. Positive focal lengths only, and no
round lens may be stronger than f = 4 mm (the objective, f = 1.5 mm).

Two parts of the real column are represented only by what they do to a round beam, as the published
descriptions of the Nion column give them. The **aberration corrector** and the quadrupole module that
couples it into the objective sit between C3 and OL1; they are not modelled, so C3 to OL1 is a plain
drift -- a parallel beam stays parallel and keeps its width. And **EL** stands in, as a round lens, for
the quadrupole-octupole coupling module that sets the EELS entrance crossover and camera length.
Element positions are fixed (edit `COLUMN` in `column.js` to change them).

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
