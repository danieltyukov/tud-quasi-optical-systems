# EE4725 - Quasi Optical Systems 2026
## Assignment 2: Parabolic Reflector Antenna

### System Parameters
- Frequency: f = 500 GHz, wavelength lambda = 0.6 mm, k0 = 2*pi/lambda
- Feed: circular aperture, diameter Df = 4*lambda = 2.4 mm (radius a = 2*lambda), y-polarized uniform electric current
- Parabolic reflector focal length: f = 5 m
- Free-space impedance: zeta = 120*pi Ohm

---

## Question 1: Feed Far Field (1 point)

**Script:** `matlab/Q1_feed_farfield.m`
**Figure:** `figures/Q1_feed_farfield.png`

The normalized far field patterns of the circular aperture feed are plotted in two principal planes:
- **E-plane** (phi = 90 deg): contains the y-axis (current polarization direction)
- **H-plane** (phi = 0 deg): perpendicular to the current direction

The far field shape is governed by the Airy diffraction pattern, arising from the 2D Fourier transform of the uniform circular aperture. The FT of a circular disk of radius a gives the characteristic 2*J1(k_perp*a)/(k_perp*a) dependence, where J1 is the first-order Bessel function and k_perp = k0*sin(theta). The computed first null at theta = 17.8 deg matches the theoretical value arcsin(1.22*lambda/Df) = arcsin(1.22/4) = 17.76 deg.

The E-plane and H-plane patterns differ due to the spectral Green's function (SGF) projection. In the H-plane (phi = 0), the SGF introduces an additional cos(theta) factor on the dominant E_theta component, causing the sidelobes to be higher than in the E-plane. Conversely, the E-plane shows deeper nulls between sidelobes. Both patterns share the same main lobe width (governed by the Airy function), but the sidelobe structure diverges at larger angles. This E/H asymmetry is a fundamental property of linearly polarized aperture antennas.

---

## Question 2: Aperture Current Distribution, f/D = 2 (2 points)

**Script:** `matlab/Q2_aperture_currents.m`
**Figures:** `figures/Q2_aperture_currents_xpol.png`, `figures/Q2_aperture_currents_ypol.png`

For f/D = 2, the reflector diameter is D = 2.5 m and the half-subtended angle is theta0 = 2*arctan(D/(4f)) = 14.25 deg. The equivalent aperture currents are computed using Love's equivalence principle:
- Magnetic surface current: Ms = -z_hat x Ea, giving Mx = Ea_y, My = -Ea_x
- Electric surface current: Js = z_hat x Ha, giving Jx = -Ea_x/zeta, Jy = -Ea_y/zeta

The aperture field is obtained by computing the feed far field at angles (theta', phi') corresponding to each aperture point via the paraboloid mapping rho = 2f*tan(theta'/2), then applying the GO amplitude taper (1+cos(theta'))/(2f).

**x-polarization figure:**
- Mx (= Ea_y, co-pol): shows a smooth, nearly circularly symmetric amplitude distribution with gradual taper from center to edge. The taper is approximately 3-4 dB at the aperture edge, caused by both the feed pattern roll-off and the GO spreading factor.
- Jx (= -Ea_x/zeta, cross-pol): displays a clear quadrupolar (four-lobed) pattern with deep nulls along both the x = 0 and y = 0 planes. The peak cross-pol level is roughly 15-20 dB below the co-pol peak. This pattern arises from the SGF coupling between the y-directed feed current and the x-component of the radiated field.

**y-polarization figure:**
- My (= -Ea_x, cross-pol): quadrupolar pattern identical in shape to Jx (both come from Ea_x).
- Jy (= -Ea_y/zeta, co-pol): smooth circular distribution identical in shape to Mx (both come from Ea_y), differing only by a factor of 1/zeta.

The co-pol components (Mx, Jy) are dominant and circularly symmetric, while cross-pol components (My, Jx) are weaker with the characteristic four-lobed pattern having nulls on the principal planes.

---

## Question 3: Reflector Far Field, f/D = 1 (2 points)

**Script:** `matlab/Q3_reflector_farfield.m`
**Figure:** `figures/Q3_reflector_farfield.png`

For f/D = 1, D = 5 m, and theta0 = 2*arctan(1/4) = 28.07 deg. The aperture currents are computed on a 301x301 grid, then the far field is obtained via numerical 2D Fourier transform using the Balanis spatial-domain formulation:
- E_theta proportional to L_phi + zeta*N_theta
- E_phi proportional to L_theta - zeta*N_phi

where N = FT(Js), L = FT(Ms). This formulation ensures proper constructive interference between the electric and magnetic current contributions in the forward direction.

The beam is extremely narrow: the theoretical Airy null is at arcsin(1.22*lambda/D) = 0.0084 deg for D = 5 m at lambda = 0.6 mm.

**Comparison with Airy pattern:**
The reflector far field pattern (both E and H planes) is compared against the Airy pattern 2*J1(k*a*sin(theta))/(k*a*sin(theta)), which corresponds to a uniformly illuminated circular aperture. The key differences are:

1. **Wider main beam**: The reflector pattern has its first null at approximately 0.01 deg, compared to 0.0084 deg for the Airy pattern. This widening is due to the amplitude taper across the aperture — the edges are less illuminated, making the effective aperture smaller.

2. **Lower sidelobes**: The first sidelobe of the reflector pattern is approximately -25 dB, significantly below the Airy pattern's first sidelobe at -17.6 dB. The amplitude taper acts as a natural window function that suppresses sidelobes.

3. **E/H plane similarity**: The E-plane and H-plane patterns are very similar, with only minor differences in the sidelobe region. This is expected because the co-pol aperture distribution is nearly circularly symmetric.

---

## Question 4: Efficiencies vs D (3 points)

**Script:** `matlab/Q4_efficiencies.m`
**Figure:** `figures/Q4_efficiencies.png`

Three efficiencies are computed as a function of reflector diameter D for 0.6 < f/D < 6 (corresponding to D from 0.83 m to 8.33 m):

**Spillover efficiency** eta_s: fraction of feed power intercepted by the reflector, computed by integrating the feed radiation intensity over the reflector's subtended solid angle (0 to theta0) and dividing by the total radiated power.

**Taper efficiency** eta_t: measures how uniformly the aperture is illuminated. Uses the Ludwig-3 co-pol definition E_co = E_theta*sin(phi) + E_phi*cos(phi) for y-polarized feed, with the formula eta_t = (16/pi)*(f/D)^2 * |integral of |E_co|*tan(theta'/2)|^2 / integral of G*sin(theta').

**Aperture efficiency** eta_ap = eta_s * eta_t.

**How f/D influences the efficiencies:**

- **Large D (small f/D, deep dish)**: The reflector subtends a large angle (theta0 up to ~45 deg at f/D = 0.6). Spillover efficiency is high (~0.96) because the large dish captures most of the feed power. However, taper efficiency is poor (~0.3) because the feed pattern has already passed through its first null at 17.8 deg and the illumination is highly non-uniform across the large aperture.

- **Small D (large f/D, shallow dish)**: The reflector subtends a small angle. Taper efficiency approaches 1.0 because only the smooth central portion of the feed pattern illuminates the dish. However, spillover efficiency drops significantly (below 0.5 at f/D > 4) because most of the feed power misses the small reflector.

- **Optimal D**: The aperture efficiency peaks at eta_ap = 0.75 at f/D = 2.28 (D = 2.19 m), where the competing effects are best balanced. At this point, eta_s = 0.87 and eta_t = 0.86. This corresponds to an edge taper of roughly -10 to -12 dB, which is the well-known criterion for maximum aperture efficiency in reflector antenna design.

---

## Question 5: Directivity and Gain vs D (2 points)

**Script:** `matlab/Q5_directivity_gain.m`
**Figure:** `figures/Q5_directivity_gain.png`

Three quantities are plotted in dBi as a function of D:

- **D_max = (pi*D/lambda)^2**: maximum directivity for a uniformly illuminated aperture
- **D = eta_t * D_max**: directivity accounting for non-uniform illumination
- **G = eta_ap * D_max = eta_s * eta_t * D_max**: gain accounting for both taper and spillover

**How f/D influences these parameters:**

- **D_max** increases monotonically with D, following the (D/lambda)^2 scaling. At D = 5 m, D_max = (pi*5/0.6e-3)^2 = 88.4 dBi. This is purely geometric — larger aperture means more collecting area.

- **Directivity D** is reduced from D_max by the taper efficiency. At large D (small f/D), the taper loss is substantial (~5-6 dB) because the feed pattern illuminates the aperture very non-uniformly. At small D (large f/D), D approaches D_max because only the smooth central part of the feed pattern is used and taper efficiency is near unity.

- **Gain G** is further reduced by spillover loss. At small D (large f/D), the gap between G and D is large (up to 5+ dB) because most feed power misses the small reflector. At large D (small f/D), G is close to D because the big dish captures nearly all the feed power (eta_s close to 1).

- The relationship **G <= D <= D_max** always holds. The gap (D_max - G) in dB equals -10*log10(eta_ap), which is minimized at the optimal f/D = 2.28 where eta_ap = 0.75 (loss of only 1.25 dB).
