# Speech Notes — Final Project 8.2

**Ka-band Satcom Dipole Array with Backing Reflector**
Daniel Tyukov | 5714699

> *Tight script: ~11 min of clean reading time, ~15 min with normal speech pace and filler words. Slow down on the figures, speed up on the bullet slides.*

---

## Slide 1 — Title  *(~15 s)*

Good morning. I'm Daniel Tyukov. I'll present my final project, number 8.2 — a Ka-band Satcom dipole array with a backing reflector.

---

## Slide 2 — Why Ka-band Satcom-on-the-move  *(~45 s)*

The application is in-flight Satcom: a flat phased array on the fuselage skin pointing at a geostationary satellite. The Ka-band transmit window is twenty-seven and a half to thirty-one gigahertz — twelve percent fractional bandwidth.

I picked a dipole over a PEC backing because it's low-profile, single-layer, and naturally resonant in this kind of bandwidth. The aircraft moves, so the array has to scan electronically, and the spec asks for a two-degree pencil beam.

---

## Slide 3 — Specifications and free design parameters  *(~35 s)*

Three hard specs: active reflection coefficient below minus ten decibels across the band, half-power beamwidth of two degrees, and no grating lobes or scan blindness.

Free parameters: dipole length and width, reflector height, lattice periods, and the feed reference impedance.

---

## Slide 4 — Approach in one slide  *(~30 s)*

I reused the periodic Method-of-Moments code from Assignments five and six. The only thing that genuinely changes is the spectral Green's function, which now has to know about the backing reflector. After that I sweep the geometry, check the scan behaviour, and size the array for the two-degree beam.

---

## Slide 5 — Theory recap: active impedance  *(~45 s)*

This is the kernel result. The active input impedance of an infinite array is a discrete Floquet sum: spectral Green's function at the Floquet wavenumber, times the squared Fourier transform of the basis function.

It splits into the fundamental mode that radiates, and a reactive higher-order tail from inter-cell coupling. Twenty-five modes per axis is converged to within one percent — same truncation I validated in Assignment five.

---

## Slide 6 — The image-theorem Green's function  *(~50 s)*

This is the only mathematical change. A horizontal current at height h above a PEC has an image at minus h with opposite sign. Adding source and image, the spectral Green's function picks up a factor of one minus the exponential of minus j twice k z h.

At broadside this becomes four sine squared of k z h, peaking at four when h is a quarter wavelength — that's where the textbook number comes from. As a free bonus, the same factor exponentially attenuates the evanescent higher-order modes, taming the inductive reactance we saw in Assignment five.

---

## Slide 7 — Lattice choice: sub-half-wave  *(~40 s)*

A first-order Floquet mode enters the visible region when sine theta-G-L equals lambda over d minus one. To keep that above unity at the highest frequency, the lattice has to be smaller than half the smallest in-band wavelength.

The smallest in-band wavelength is just under ten millimetres, so I picked four point eight. Lambda over d stays above two everywhere in the band — no grating lobes at any scan angle, and as a bonus, no Wood-Rayleigh scan-blindness either.

---

## Slide 8 — Q1 design sweep at broadside  *(~60 s)*

I swept dipole length, reflector height, and feed impedance, and minimised the worst-case in-band reflection.

Top-left: input impedance. The real part crosses one hundred ohms in the centre of the band. Top-right: the active reflection coefficient is minus seventeen and a half decibels in the worst part of the band — seven and a half decibels of margin below spec. The red curve shows what happens without the reflector — barely minus three decibels, well out of spec.

Bottom-left contour shows the design region; my point sits comfortably inside it. Bottom-right confirms one hundred ohms is the optimum feed reference.

---

## Slide 9 — The reflector is what makes the bandwidth  *(~45 s)*

Two design choices worth flagging. First: I'm using h equals one fifth of a wavelength, not a quarter. Quarter wavelength makes the resonance sharp; a shorter reflector flattens it and broadens the band, which I need for twelve percent bandwidth.

Second: the dipole is forced up against the cell boundary at four point five millimetres. A half-wave dipole would be five point one — too long for the cell. So the dipole resonance lands just above the centre frequency, and Re-Z stays in the matching window across the band.

---

## Slide 10 — Q2 grating-lobe diagrams  *(~40 s)*

Visual confirmation. Each panel is a grating-lobe diagram at one in-band frequency. Blue disc is the visible region; grey circles are the higher-order Floquet circles.

Even at thirty-one gigahertz, where lambda over d is at its smallest, the grey circles still sit outside the visible region. The dashed lines are the E and H scan trajectories — they stay entirely inside. So: no grating lobes anywhere in the band.

---

## Slide 11 — Q2 scan study  *(~55 s)*

Active reflection coefficient versus scan angle, in all three principal planes, at the band edges and centre.

The dashed red lines mark where the worst case across frequencies crosses minus ten decibels. E-plane gives twenty-eight degrees, H-plane thirty-seven, D-plane forty-three. The system limit is the worst plane and worst frequency — twenty-eight degrees in the E-plane at the high band edge.

The curves are smooth — no Wood-anomaly spikes — exactly because the lattice is sub-half-wave. The roll-off is just smooth scan loss.

---

## Slide 12 — Why the E-plane is the bottleneck  *(~45 s)*

Worth a quick comment, because students usually expect the H-plane to be worse.

H-plane: scan moves along k-y, basis transform is unchanged, only one over cos theta in the wave impedance — so resistance grows with scan.

E-plane: scan moves along k-x, which directly attacks the basis Fourier transform, and there's a cos-theta factor on top. Re-Z drops below the matching window first. That's why the E-plane limit is twenty-eight degrees — about ten less than the H-plane.

---

## Slide 13 — Q3 array sizing  *(~50 s)*

Two-degree half-power beamwidth. From the closed-form expression — zero point eight eight six lambda over array length — the worst case is the lowest frequency, where lambda is largest. The table gives fifty-eight elements at twenty-seven and a half gigahertz, fifty-five at the centre, fifty-two at thirty-one.

The figure is the numerical check: simulation overlaps the closed form perfectly. Fifty-five per side gives one point nine seven degrees at the centre — just under spec. So my recommendation is fifty-five by fifty-five, three thousand and twenty-five elements. Fifty-eight by fifty-eight if you want strict worst-case across the whole band.

---

## Slide 14 — Q3 beamwidth across the band  *(~30 s)*

Just for context. Beamwidth scales as one over frequency, so the fifty-five-by-fifty-five design narrows to one point eight five degrees at the high band edge and widens to two point zero nine at the low edge — a four percent overshoot, which I consider a non-issue.

---

## Slide 15 — Conclusions  *(~45 s)*

To wrap up. One mathematical change — the image-theorem Green's function — drives the entire result. The chosen unit cell hits minus seventeen and a half decibels at broadside, seven and a half below spec. The sub-half-wave lattice rules out grating lobes and scan blindness. The system scan limit is plus or minus twenty-eight degrees, set by the E-plane at the high band edge. And fifty-five by fifty-five elements deliver the two-degree pencil beam.

Thank you — happy to take questions.
