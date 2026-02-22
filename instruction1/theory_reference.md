# EE4580 Quasi Optical Systems - Topic 1: Free Space Spectral Green's Function

## Course Info
- Course: EE4580 / ET4580, TU Delft
- Teachers: A. Neto, D. Cavallo, N. Llombart (THz Sensing Group)
- Week 1 topic: Spectral Green's Function (A. Neto)
- Structure: 2h theory + 2h MATLAB instruction per week
- Deadline: Sunday 22nd at 23:59

## Learning Objectives
1. Free Space Green's function
2. Spatial Fourier Transforms
3. Asymptotic evaluation of integrals via saddle point

---

## Part I: Green's Functions in Spatial Domain

### What is a Green's Function?
- G(r - r') is the field G at observation point r, radiated by an elementary equivalent current source delta(r - r') * p_hat at source point r'.
- Electric field from electric currents: e(r) = integral over V of G^ej(r - r') * j(r') dr'
- This is a convolution integral between the GF and equivalent currents.

### Free Space Dyadic GF (spatial domain)
For a z-oriented point source j(r) = delta(x)delta(y)delta(z) z_hat:
- Magnetic vector potential: A(r,j) = mu * e^{-jk|r-r'|} / (4pi|r-r'|) * z_hat
- E field from potential: e(r) = -j*omega*A - j/(omega*epsilon*mu) * grad(div(A))
- H field from potential: h(r) = (1/mu) * curl(A)

Far field of z-oriented Hertzian dipole:
- h_far = jk * sin(theta) * e^{-jkr} / (4*pi*r) * phi_hat
- e_far = jk*zeta * sin(theta) * e^{-jkr} / (4*pi*r) * theta_hat

General dyadic GF (for arbitrary dipole orientation p_hat):
```
G^ej_fs(r-r') = -j*(zeta/k) * [differential_matrix] * e^{-jk|r-r'|} / (4*pi*|r-r'|)
```
where the differential matrix is:
```
| k^2 + d^2/dx^2    d^2/dxdy       d^2/dxdz   |
| d^2/dydx          k^2 + d^2/dy^2 d^2/dydz   |
| d^2/dzdx          d^2/dzdy       k^2+d^2/dz^2|
```
Problem: derivatives still need to be performed (not fully analytical).

---

## Part II: Spectral Representation

### Fourier Transform Convention
```
G(kx) = integral g(x) * e^{j*kx*x} dx        (FT)
g(x)  = 1/(2pi) * integral G(kx) * e^{-j*kx*x} dkx  (IFT)
```

### Dirac Delta in Spectral Domain
- FT{delta(x-x')} = e^{j*kx*x'}
- IFT{e^{j*kx*x'}} = delta(x-x')
- 3D: delta(r-r') = 1/(2pi)^3 * triple_integral e^{-jkx(x-x')} e^{-jky(y-y')} e^{-jkz(z-z')} dkx dky dkz

### Helmholtz Equation and Scalar Potential
```
nabla^2 * psi + k^2 * psi = -delta(r - r')
```
Solution (with radiation condition): psi = e^{-jk|r-r'|} / (4pi|r-r'|)

### Spectral Representation of Scalar GF (3D triple integral)
```
e^{-jk|r-r'|} / (4pi|r-r'|) = 1/(2pi)^3 * triple_integral
    -e^{-jkx(x-x')} e^{-jky(y-y')} e^{-jkz(z-z')} / (k^2 - kx^2 - ky^2 - kz^2) dkx dky dkz
```

Derivation steps:
1. Express psi via its anti-Fourier transform: psi = 1/(2pi)^3 * integral Psi(k_vec, r') e^{-jkx*x} e^{-jky*y} e^{-jkz*z} dkx dky dkz
2. Substitute into Helmholtz equation
3. In spectrum: nabla^2 -> -(kx^2 + ky^2 + kz^2)
4. Equate spectra (valid for every x,y,z):
   Psi(kx,ky,kz,r') = -e^{jkx*x'} e^{jky*y'} e^{jkz*z'} / (k^2 - kx^2 - ky^2 - kz^2)

### Spectral Derivatives (key simplification!)
In spectral domain, spatial derivatives become algebraic:
- d^2/dx^2 -> -kx^2
- d^2/dxdy -> -kx*ky
- d^2/dxdz -> -kx*kz
- etc.

### Spectral GF (3D triple integral form)
```
G^ej_fs(r-r') = j*(zeta/k) * 1/(2pi)^3 * triple_integral
    [spectral_matrix] * e^{-jkx(x-x')} e^{-jky(y-y')} e^{-jkz(z-z')} / (k^2 - kx^2 - ky^2 - kz^2) dkx dky dkz
```
where spectral_matrix:
```
| k^2 - kx^2    -kx*ky      -kx*kz  |
| -ky*kx         k^2 - ky^2 -ky*kz  |
| -kz*kx         -kz*ky      k^2-kz^2|
```

---

## Part III: 2D Integral Representation (key result)

### Spectral Identity (reducing 3D to 2D integral)
For any (z - z') != 0:
```
e^{-jk|r-r'|} / (4pi|r-r'|) = -j/(8pi^2) * double_integral
    e^{-jkx(x-x')} e^{-jky(y-y')} e^{-j*sqrt(k^2-kx^2-ky^2)*|z-z'|} / sqrt(k^2 - kx^2 - ky^2) dkx dky
```

Define: kz = sqrt(k^2 - kx^2 - ky^2)

### Derivation via Complex Contour Integration
The kz integral is evaluated using residue calculus:
```
integral_{-inf}^{inf} e^{-jkz(z-z')} / ((kp - kz)(kp + kz)) dkz = j*pi * e^{-jkp|z-z'|} / kp
```
where kp^2 = k^2 - kx^2 - ky^2

Key steps:
- Two poles at kz = +/- kp
- For (z-z') > 0: close contour in lower half-plane (Im(kz) < 0), enclose pole at +kp
- For (z-z') < 0: close contour in upper half-plane (Im(kz) > 0), enclose pole at -kp
- Both cases give same result with |z-z'|

### Complex Plane Topology
- k = kr - j*ki (lossy medium, ki >= 0)
- kp = kpr - j*kpi with kpr > 0, kpi > 0
- Poles at +kp (4th quadrant) and -kp (2nd quadrant) in kz plane

### Generalized Spectral Identity
For analytical functions C(kz):
```
integral C(kz) * e^{-jkz(z-z')} / ((kp-kz)(kp+kz)) dkz = j*pi * e^{-jkp|z-z'|} / kp * C(+/-kp)
```
where +kp for (z-z') > 0, -kp for (z-z') < 0

### 2D Spectral GF (FINAL KEY RESULT)

**For electric field due to electric currents:**
```
G^ej_fs(r,r') = -(zeta/k) * 1/(8*pi^2) * double_integral
    [D_ej matrix] * e^{-jkx(x-x')} e^{-jky(y-y')} e^{-jkz|z-z'|} / kz dkx dky
```
where kz = sqrt(k^2 - kx^2 - ky^2), and the +/- kz in the matrix depends on sign of (z-z'):
```
D_ej = | k^2 - kx^2     -kx*ky         -kx*(+/-kz)    |
       | -ky*kx          k^2 - ky^2     -ky*(+/-kz)    |
       | -(+/-kz)*kx     -(+/-kz)*ky    k^2 - kz^2     |
```

**For electric field due to magnetic currents:**
```
G^em_fs(r,r') = -j/(8*pi^2) * double_integral
    [D_em matrix] * e^{-jkx(x-x')} e^{-jky(y-y')} e^{-jkz|z-z'|} / kz dkx dky
```
```
D_em = | 0           +/-j*kz    -j*ky |
       | -/+j*kz     0           j*kx |
       | j*ky        -j*kx       0    |
```

**For magnetic field due to magnetic currents:**
Same structure as G^ej but with 1/zeta instead of zeta.

**For magnetic field due to electric currents (G^hj):**
```
G^hj_fs(r,r') = j/(8*pi^2) * double_integral
    [D_hj matrix] * e^{-jkx(x-x')} e^{-jky(y-y')} e^{-jkz|z-z'|} / kz dkx dky
```

### Special case z = z'
The zz-component of the GF contains a delta function term:
```
(k^2 - kz^2 - kp^2 + kp^2)/(k^2 - kp^2 - kz^2) integrand over kz
= delta(z-z') + j/2 * (k^2-kz^2)/sqrt(k^2-kp^2) * e^{-j*sqrt(k^2-kp^2)*|z-z'|}
```

### Key Physical Interpretation: Plane Wave Expansion
- The spectral GF is a plane wave representation of EM fields
- Each (kx, ky) pair defines a plane wave direction
- kx^2 + ky^2 + kz^2 = k^2 (dispersion relation)
- Space dependence is only in the exponential
- For kx^2+ky^2 < k^2: propagating waves (real kz)
- For kx^2+ky^2 > k^2: evanescent waves (imaginary kz)

### Convergence Notes
- Diagonal terms converge slowly (sum of potential term + derivative terms)
- For r = r', GF diverges (not an observable quantity alone, only meaningful in convolution)

---

## Part IV: Far Field and Stationary Phase

### Radiation from Current Distributions
For planar source at z = z':
```
e(x,y,z) = 1/(4pi^2) * double_integral D^fc(kx,ky) * J_tilde(kx,ky) * e^{-jkx*x} e^{-jky*y} * e^{-jkz|z-z'|}/kz dkx dky
```
where:
- D^fc is one of the four 3x3 spectral dyads
- J_tilde(kx,ky) = FT of the current distribution = integral j(x',y') e^{jkx*x'+jky*y'} dx'dy'

### Far Field via Stationary Phase Point
At the stationary phase point (kxs, kys, kzs), the slowly-varying part of the integrand can be extracted:
```
f^far(r) = j * D^fc(kxs, kys, z, z') * C_tilde(kxs, kys) * e^{-jkr} / (2*pi*r)
```

The stationary phase point corresponds to the direct ray from source to observation:
- kxs = k*sin(theta)*cos(phi)
- kys = k*sin(theta)*sin(phi)
- kzs = k*cos(theta)

### Stationary Phase Method (Appendix)
For integrals of the form I = integral f(x) * e^{-j*Omega*q(x)} dx:
1. Find stationary phase point x0 where q'(x0) = 0
2. Taylor expand q(x) around x0: q(x) ~ q(x0) + q''(x0)/2 * (x-x0)^2
3. Result: I ~ f(x0) * e^{-j*Omega*q(x0)} * sqrt(pi / (j*Omega*q''(x0)/2))

Key identity: integral e^{-j*Omega*x^2} dx = sqrt(pi/|Omega|) * e^{∓j*pi/4}