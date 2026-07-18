# Assignment 1: Farfield Evaluation Using SGF

## Course & Logistics
- **Course:** EE4580 Quasi Optical Systems, TU Delft
- **Topic:** MATLAB Instruction - Farfield evaluation using SGF
- **Deadline:** Sunday 22nd at 23:59
- **Submission:** PDF report (exported from .mlx) + all MATLAB code (.mlx + .m files)
- **Important:** Code must run without errors on a clean MATLAB session. Restart MATLAB and execute from beginning before submitting.

---

## Theory: Farfield Evaluated Using SGF

### Electric field from current (spectral domain convolution)
```
e(x,y,z) = 1/(4*pi^2) * double_integral G_ej(kx,ky) * J_tilde(kx,ky) * e^{-j*kx*x} * e^{-j*ky*y} * e^{-j*kz*|z-z'|} dkx dky
```

### Far field via stationary phase (asymptotic evaluation)
```
E^far(r, theta, phi) = j * kzs * G_ej(kxs, kys, z, z') * J_tilde(kxs, kys) * e^{-j*k*r} / (2*pi*r)
```
where:
- (x', y', z') = source position, (x, y, z) = observation position
- kxs = k * sin(theta) * cos(phi)
- kys = k * sin(theta) * sin(phi)
- kzs = k * cos(theta)
- J_tilde(kx, ky) = FT of the current distribution
- G_ej = dyadic SGF for electric currents

### Dyadic SGF for electric currents
```
G_ej(kx, ky) = -(zeta) / (2*k*kz) * [ k^2-kx^2    -kx*ky      -kx*(+/-kz) ]
                                       [ -kx*ky      k^2-ky^2    -ky*(+/-kz) ]
                                       [ -kx*(+/-kz) -ky*(+/-kz)  k^2-kz^2   ]
```
where +kz for z > z', -kz for z < z'.

### kz definition (unique)
```
kz = -j * sqrt( -(k^2 - kx^2 - ky^2) )
```
This ensures proper behavior: real kz for propagating waves (kx^2+ky^2 < k^2), imaginary kz with correct sign for evanescent waves.

### Directivity
```
D(theta, phi) = U(theta, phi) / (Prad / (4*pi))
```
where radiation intensity:
```
U(theta, phi) = |E^far(r, theta, phi)|^2 / (2*zeta) * r^2
```
and total radiated power:
```
Prad = integral_0^{2*pi} integral_0^{pi} U(theta, phi) * sin(theta) d(theta) d(phi)
```

---

## MATLAB Implementation Steps (Reusable Routines)

### Step 1: Calculate the dyadic SGF
```matlab
[ej_SGF] = EJ_SGF(er, k, kx, ky)
% Returns 3x3 matrix (or set of components) of the spectral dyadic GF
% er = relative permittivity (1 for free space)
% k = wavenumber
% kx, ky = spectral variables
```

### Step 2: FT of the current distribution
```matlab
Jx = FTCurrent(k, er, kx, ky, l, w)
% k = wavenumber (k0, free-space, per professor's clarification)
% er = relative permittivity
% kx, ky = spectral variables
% l = dipole length, w = dipole width
```

### Step 3: Evaluate the far field
```matlab
[Eth, Eph] = farfield(k, R_FF, TH, PH, kz, Gxx, Gyx, Gzx, Jx)
% k = wavenumber
% R_FF = far field distance
% TH, PH = theta and phi angles
% kz = z-component of wavevector at stationary phase point
% Gxx, Gyx, Gzx = relevant SGF components
% Jx = FT of x-directed current
```

### Step 4: Estimate the directivity
```matlab
[Dir, Prad] = Directivity(E_tot, Theta, dth, dph, er, r)
% E_tot = total electric field magnitude
% Theta = theta angle array
% dth, dph = angular step sizes
% er = relative permittivity
% r = observation distance
```

---

## Q1 (3 points): Elementary Electric Source

### Setup
- x-oriented elementary electric source at origin
- j(x,y,z) = delta(x)*delta(y)*delta(z) * x_hat
- FT of this current: J_tilde(kx, ky) = x_hat (just a unit vector)
- Frequency: f = 30 GHz
- Spectral range: kx in [0, 3*k0], ky = 0

### Q1.1: Plot x-, y-, z-components of SGF
- Plot the x, y, and z components of the SGF given the x-oriented current
- Both **real and imaginary** parts should be plotted
- Since current is x-directed, the relevant SGF column is the first column: Gxx, Gyx, Gzx

### Q1.2: Imaginary region of x-component
- Identify in which region of the spectrum the x-component is imaginary
- Justify with the plot and the SGF expression
- (Answer relates to: when kx^2 + ky^2 > k^2, kz becomes imaginary, changing which parts of Gxx are real vs imaginary)

---

## Q2 (4 points): Farfield of a Dipole

### Dipole geometry
- x-directed dipole centered at origin
- Oriented along x-axis with length L and width W along y
- Radiating in free space (upper hemisphere, z > 0)

### Current distribution model (open-circuited transmission line approximation)
The spatial current on the dipole:
```
Jx(x, y) = l(x) * t(y) * x_hat
```
Its Fourier transform:
```
J_tilde_x^FT(kx, ky) = L(kx) * T(ky)
```

#### Transverse distribution (y-direction):
```
t(y) = (1/w) * rect(y, w)

T(ky) = sinc(ky * w / 2)
```
Note: sinc here likely means sin(x)/x (unnormalized sinc).

#### Longitudinal distribution (x-direction):
```
l(x) = sin(keq * (l/2 - |x|)) / sin(keq * l/2)

L(kx) = 2*keq * (cos(kx*l/2) - cos(keq*l/2)) / ((keq^2 - kx^2) * sin(keq * l/2))
```

#### keq definition
**Professor's clarification:** Ignore the original keq definition in this exercise. Use just the free-space wavenumber:
```
keq = k0    (NOT sqrt((kx^2 + ky^2)/2) as shown in the original)
```

**Professor's clarification:** Use the analytical formulas for L(kx) and T(ky) directly (do not numerically compute the FT).

### Q2.1: 1D plot of farfield in main planes
- Dipole dimensions: L = lambda0/2, w = lambda0/40 at 30 GHz
- Calculate far field (total field) for phi = 0deg, 45deg, 90deg
- Plot as function of theta

### Q2.2: Farfield in upper medium (z>0) using UV representation
- UV plane: U = sin(theta)*cos(phi), V = sin(theta)*sin(phi)
- Use 'surface' (top view), 'pcolor', or 'imagesc' plot
- Dynamic range of 10 dB

### Q2.3: Directivity at broadside in dB vs frequency
- Fixed dipole dimensions: L = 5 mm, W = 0.25 mm
- Calculate far field as function of frequency
- Evaluate directivity for each frequency
- Plot directivity at broadside (theta=0) in dB vs frequency

---

## Q3 (3 points): Dipole with a Backing Reflector (PEC ground plane)

### Setup
- Same dipole as Q2: L = lambda0/2, w = lambda0/40 at 30 GHz
- Dipole at height h above an infinite PEC plane
- Use image theorem to replace PEC with image sources in free space

### Image theorem for PEC ground plane
For an electric conductor (n_hat x E = 0):
- **Electric current J:** image is J' with **same vertical component, opposite horizontal component**
  - Real source: J pointing right (x-directed) -> Image: J' pointing LEFT (-x-directed)
  - The image is at z = -h (mirrored below PEC at z = 0)
- **Magnetic current M:** image is M' with **opposite vertical component, same horizontal component**

So for x-directed dipole at height h:
- Real source at z = +h: current +x_hat
- Image source at z = -h: current -x_hat (opposite sign for tangential electric current)

The total far field = contribution from real source + contribution from image source, both computed using free-space SGF.

### Q3.1: 1D farfield along main planes, compare to Q2
- h = 7.5 mm
- Plot farfield in main planes (phi = 0deg, 45deg, 90deg)
- Compare with Q2 results (free-space dipole)

Before computing numerically, write down the analytical expression of the current spectrum.

**Steps:**
1. Calculate the dyadic SGF
2. Calculate the FT of the current distribution (now includes both real + image source)
3. Evaluate the far field (total field)

### Q3.3: Directivity at broadside (linear scale) vs distance h from ground plane
- Plot from h = lambda0/10 to h = 2*lambda0
- x-axis in terms of wavelength (h/lambda0)

### Additional questions for Q3.3:
1. **Integration domain for radiated power:** In which domain do you integrate? (Hint: since image theorem is applied, integrate over upper hemisphere only, i.e., 0 to pi/2 for theta, since the ground plane means no radiation below)
2. **Null at broadside:** For which values of h do you expect a null at broadside? (Think about kz expression - when the real and image contributions cancel at theta=0)
3. **Power radiated for h = 0:** What is Prad when h = 0? Use the spectrum and/or image to explain. Verify numerically. (When h=0, the image current -x_hat coincides with the real current +x_hat, so they cancel and Prad = 0)

### Important note from assignment:
**"Please explain all the figures and results!"** - Every plot needs accompanying discussion/explanation.

---

## Professor's Suggestions & Clarifications (Posted 16 Feb 2026)

1. **keq in Q2:** Ignore the original definition. Use just the free-space wavenumber k0.
2. **FT of current on dipole:** Use the analytical formulas provided (L(kx) and T(ky)) directly. Do not compute the FT numerically.
3. **Far-field evaluation:** Use the asymptotic (stationary phase) evaluation:
   ```
   E^far(r, theta, phi) = j * kzs * G_ej(kxs, kys, z, z') * J_tilde(kxs, kys) * e^{-j*k*r} / (2*pi*r)
   ```

---

## Key Constants & Parameters Summary

| Parameter | Value | Notes |
|-----------|-------|-------|
| f | 30 GHz | Operating frequency |
| lambda0 | c/f = 10 mm | Free-space wavelength at 30 GHz |
| k0 | 2*pi/lambda0 | Free-space wavenumber |
| zeta | 120*pi ohms (~377 ohms) | Free-space impedance |
| c | 3e8 m/s | Speed of light |
| er | 1 | Free space relative permittivity |

### Q1 parameters
| Parameter | Value |
|-----------|-------|
| kx range | [0, 3*k0] |
| ky | 0 |
| Current | delta source, x-oriented |

### Q2 parameters
| Parameter | Q2.1/Q2.2 | Q2.3 |
|-----------|-----------|------|
| L | lambda0/2 = 5 mm | 5 mm (fixed) |
| W | lambda0/40 = 0.25 mm | 0.25 mm (fixed) |
| phi cuts | 0, 45, 90 deg | broadside only |
| keq | k0 (clarified) | k0 (clarified) |

### Q3 parameters
| Parameter | Q3.1 | Q3.3 |
|-----------|------|------|
| L | lambda0/2 = 5 mm | lambda0/2 = 5 mm |
| W | lambda0/40 = 0.25 mm | lambda0/40 = 0.25 mm |
| h | 7.5 mm (= 3*lambda0/4) | lambda0/10 to 2*lambda0 |
| PEC plane | at z = 0 | at z = 0 |

---

## File Structure (Recommended)
```
instruction1/
  EJ_SGF.m          - Function: spectral dyadic GF
  FTCurrent.m        - Function: FT of dipole current distribution
  farfield.m         - Function: far field computation
  Directivity.m      - Function: directivity calculation
  sgf_homework1.mlx  - Live script: main assignment (all Q1-Q3, plots, explanations)
  assignment1_details.md  - This reference file
```
