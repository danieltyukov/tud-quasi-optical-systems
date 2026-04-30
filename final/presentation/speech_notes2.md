# Speech Notes — Final Project 8.2

**Ka-band Satcom Dipole Array with Backing Reflector**
Daniel Tyukov | 5714699

> *Trimmed for ~9 min clean reading, so with normal pace and filler words it lands comfortably under 15. Slow on figure slides (8, 11, 13). Fast on bullet slides (3, 4, 7, 14).*

---

## Slide 1 — Title  *(~10 s)*

Good morning. I'm Daniel Tyukov, and I'll present my final project — number 8.2 — a Ka-band Satcom dipole array with a backing reflector.

---

## Slide 2 — Why Ka-band Satcom-on-the-move  *(~35 s)*

The setting on the right is in-flight Satcom — a flat phased array bonded to the aircraft skin, talking to a geostationary satellite while the plane is in motion. The transmit window runs from twenty-seven and a half to thirty-one gigahertz, so roughly twelve percent fractional bandwidth.

I'm using printed dipoles over a PEC backing — single-layer, low-profile, and broadband enough for that window. Since the aircraft moves, the array has to scan electronically, and the spec at the bottom asks for a two-degree pencil beam with the active reflection coefficient below minus ten decibels across the band.

---

## Slide 3 — Specifications and free design parameters  *(~25 s)*

On the left are the three hard requirements: active reflection below minus ten decibels, half-power beamwidth of two degrees, and no grating lobes or scan blindness anywhere we plan to scan. On the right are the knobs I'm allowed to turn — dipole length, width, reflector height, the two lattice periods, and the feed reference impedance. The little side-view shows how those parameters fit together.

---

## Slide 4 — Approach in one slide  *(~20 s)*

Five steps. I reuse the periodic Method-of-Moments code from Assignments five and six, replace the free-space Green's function with one that includes the PEC image, sweep the design space at broadside, verify scanning behaviour, then size the array for the beamwidth target.

---

## Slide 5 — Theory recap: active impedance  *(~40 s)*

The boxed equation at the top is the kernel result. Galerkin plus Floquet reduces the integral equation to a discrete sum: spectral Green's function at the Floquet wavenumber, times the squared Fourier transform of the basis function, summed over all modes.

On the right is the decomposition I'll use throughout. The fundamental mode Z-zero-zero is the part that radiates into the main beam. Z-h-o is the sum over higher-order evanescent modes — that's the reactive near-field coupling between cells. Twenty-five modes per axis is converged to within one percent.

---

## Slide 6 — The one new ingredient: image-theorem Green's function  *(~45 s)*

This is the only mathematical change. A horizontal current at height h above a PEC has an image at minus h with opposite sign. Adding source and image, the spectral Green's function picks up a factor of one minus the exponential of minus j twice k z h.

At broadside this becomes four sine squared of k z h, peaking at four when h is a quarter wavelength — that's where the textbook number comes from. As a free bonus, the same factor exponentially attenuates the evanescent higher-order modes, taming the inductive reactance we saw in Assignment five.


---

## Slide 7 — Lattice choice: sub-half-wave to defeat grating lobes  *(~30 s)*

A first-order Floquet mode enters the visible region when sine theta-G-L equals lambda over d minus one. To keep that above unity at the highest frequency, the lattice has to be smaller than half the smallest in-band wavelength.

The smallest in-band wavelength is just under ten millimetres, so I picked four point eight. Lambda over d stays above two everywhere in the band — no grating lobes at any scan angle, and as a bonus, no Wood-Rayleigh scan-blindness either.


---

## Slide 8 — Q1 design sweep at broadside  *(~65 s)*

I swept dipole length, reflector height, and feed impedance, and minimised the worst-case in-band reflection.

Top-left: input impedance. The real part crosses one hundred ohms in the centre of the band. Top-right: the active reflection coefficient is minus seventeen and a half decibels in the worst part of the band — seven and a half decibels of margin below spec. The red curve shows what happens without the reflector — barely minus three decibels, well out of spec.

Bottom-left contour shows the design region; my point sits comfortably inside it. Bottom-right confirms one hundred ohms is the optimum feed reference.

---

## Slide 9 — The reflector is what makes the bandwidth  *(~35 s)*

The left box collects the final unit-cell numbers, and at the bottom you can see the worst-case reflection — minus seventeen and a half decibels with the reflector, only minus three without it.

Two design choices on the right worth flagging. First, h is one over five point five of a wavelength — a shorter reflector than a quarter wavelength flattens the resonance and broadens the band, which is what I need for twelve percent bandwidth. Second, the dipole length is ninety-five percent of the cell — pushing it as long as the geometry allows places the dipole resonance just above centre frequency, keeping the resistance in the matching window across the band.

---

## Slide 10 — No grating lobes anywhere in the band  *(~30 s)*

Visual confirmation of the lattice choice. Each panel is a grating-lobe diagram at one frequency in the band. The blue disc is the visible region, and the grey circles are the higher-order Floquet centres.

You can see that even at thirty-one gigahertz, in the right-hand panel where lambda over d is at its smallest, the visible region still doesn't touch any of the grey circles. So no grating lobes anywhere in the band.

---

## Slide 11 — Q2 scan study — three planes, three frequencies  *(~50 s)*

Active reflection coefficient versus scan angle, in all three principal planes, at the band edges and centre.

The dashed red lines mark where the worst case across frequencies crosses minus ten decibels. E-plane gives twenty-eight degrees, H-plane thirty-seven, D-plane forty-three. The system limit is the worst plane and worst frequency — twenty-eight degrees in the E-plane at the high band edge.

The curves are smooth — no Wood-anomaly spikes — exactly because the lattice is sub-half-wave. The roll-off is just smooth scan loss.

---

## Slide 12 — Why is the E-plane the bottleneck?  *(~40 s)*

The boxed equation at the top is the fundamental Floquet term of the active impedance on the principal scan trajectory. The two columns specialise it to the H-plane and E-plane.

In the H-plane, kₓ-zero stays at zero, so the basis Fourier transform doesn't see the scan at all. The only effect is one over cosine theta in the wave impedance — the resistance just grows smoothly with scan, easy to keep matched.

In the E-plane, the scan moves kₓ-zero directly into the basis Fourier transform, which rolls off, and there's a cosine theta on top. The resistance drops fast, and the match hits minus ten decibels first. That's the mechanism behind the twenty-eight-degree limit.

---

## Slide 13 — Q3 array sizing for two-degree HPBW  *(~50 s)*

Two-degree half-power beamwidth. From the closed-form expression — zero point eight eight six lambda over array length — the worst case is the lowest frequency, where lambda is largest. The table gives fifty-eight elements at twenty-seven and a half gigahertz, fifty-five at the centre, fifty-two at thirty-one.

The figure is the numerical check: simulation overlaps the closed form perfectly. Fifty-five per side gives one point nine seven degrees at the centre — just under spec. So my recommendation is fifty-five by fifty-five, three thousand and twenty-five elements. Fifty-eight by fifty-eight if you want strict worst-case across the whole band.

---

## Slide 14 — Beamwidth across the band  *(~20 s)*

Just for context. Beamwidth scales as one over frequency, so the fifty-five-by-fifty-five design narrows to one point eight five degrees at the high band edge and widens to two point zero nine at the low edge — a four percent overshoot, which I consider a non-issue.

---

## Slide 15 — Conclusions  *(~30 s)*

To wrap up. One mathematical change — the image-theorem Green's function — drove the entire result. The unit cell hits minus seventeen and a half decibels at broadside, seven and a half below spec. The sub-half-wave lattice rules out grating lobes and scan blindness. The system scan limit is plus or minus twenty-eight degrees, set by the E-plane at the high band edge. And fifty-five by fifty-five elements deliver the two-degree pencil beam.

Thank you — happy to take questions.
