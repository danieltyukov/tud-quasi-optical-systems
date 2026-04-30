# Final Project 8.2 — Study Guide

**Ka-band Satcom Dipole Array with Backing Reflector**
Daniel Tyukov | 5714699

> Comprehensive theory + design rationale + curveball-question prep for the defence.

---

## 1. What you actually did, in human language

You designed an **antenna array for an aircraft talking to a satellite**. The plane moves, so the beam has to be steered electronically. The frequency band is **27.5–31 GHz** (Ka-band, 12% bandwidth). The radiator is a **printed dipole with a metal sheet behind it** — the metal sheet ("backing reflector") makes the antenna radiate forward instead of equally in both directions.

Your job had three parts:

1. **Q1 — broadside design.** Pick the geometry (dipole length, reflector height, lattice spacing, feed impedance) so that when the array points straight up, the antenna is well-matched ($|\Gamma| < -10$ dB) across the whole band.
2. **Q2 — scanning behaviour.** Once the geometry is fixed, find out how far the beam can be steered before either the match collapses below $-10$ dB, or grating lobes show up, or scan blindness happens.
3. **Q3 — array sizing.** Given a 2° pencil-beam spec, how many elements do we need on each side?

Your computational engine is **periodic Method of Moments**. You analyse one unit cell, assume infinite periodicity, and then — for far-field patterns only — multiply by an array factor for a finite $N \times N$ truncation.

---

## 2. The theory backbone — concept by concept

### 2.1 Green's functions (the foundation)

A **Green's function** $\bar{\bar{g}}(\vec{r},\vec{r}')$ is the field at observation point $\vec{r}$ produced by a unit point current at source point $\vec{r}'$. Once you know $g$, *any* current distribution $\vec{j}(\vec{r})$ gives a field via convolution:

$$\vec{e}(\vec{r}) = \iiint_V \bar{\bar{g}}^{ej}(\vec{r}-\vec{r}') \, \vec{j}(\vec{r}') \, d\vec{r}'$$

So the GF is the antenna problem's "unit response."

**Free-space spectral GF** (dyadic, electric current → electric field, $xx$ component, evaluated in the $(k_x, k_y)$ plane):

$$G_{xx}^{ej}(k_x, k_y) = -\frac{\zeta}{2 k_0} \cdot \frac{k_0^2 - k_x^2}{k_z}$$

with $k_z = \sqrt{k_0^2 - k_x^2 - k_y^2}$. When $k_z$ is real → propagating plane wave. When imaginary → evanescent (decays in $z$).

**Why we use the spectral GF, not the spatial one.** In free space the spatial GF $e^{-jkr}/(4\pi r)$ is fine. But the moment you add a ground plane, dielectric slab, multilayer stack — the spatial GF becomes intractable, while the spectral GF just gets new multiplicative factors. Spectral GF = plane-wave decomposition. Reflections, transmissions, image theorem all work plane-wave by plane-wave.

---

### 2.2 The image theorem — your one piece of new physics

The whole final project hinges on this one trick. A horizontal electric current at height $h$ above a perfect electric conductor (PEC) is equivalent — for the field above the PEC — to that current plus a mirror image current at depth $-h$ with **opposite sign** (parallel components flip; perpendicular components don't).

When you add the original and image plane waves in the spectral domain, you get a multiplier on the free-space GF:

$$\boxed{\,G_{\mathrm{eff}}(k_x, k_y) = G_{\mathrm{FS}}(k_x, k_y) \cdot \left(1 - e^{-j\,2\,k_z\,h}\right)\,}$$

**At broadside** ($k_x = k_y = 0$, so $k_z = k_0$), the magnitude squared of that factor becomes:

$$\left|1 - e^{-j 2 k_0 h}\right|^2 = 4 \sin^2(k_0 h)$$

That's $4$ when $k_0 h = \pi/2$, i.e. $h = \lambda/4$. That's the textbook quarter-wave reflector rule.

**Bonus effect.** For higher-order Floquet modes, $k_z$ is *imaginary* (those modes are evanescent). Then $e^{-j\,2\,k_z\,h}$ becomes a real decaying exponential — so the factor $\bigl(1 - e^{-2|k_z|h}\bigr)$ grows smaller and **exponentially kills the inductive reactance** that bedevilled the unbacked array. You bought yourself two improvements with one ground plane.

---

### 2.3 Periodic Green's function & Floquet modes

**Floquet's theorem**: in a periodic structure with progressive phase shift, the field has the form $(\text{periodic function}) \times e^{-j(k_{x0} x + k_{y0} y)}$. Combined with the **Poisson summation formula**, this turns the continuous spectral integrals into discrete sums:

$$\boxed{\,\tilde{\bar{g}}_{\infty}^{ej}(x,y,z) = \frac{1}{d_x d_y} \sum_{m_x = -\infty}^{\infty} \sum_{m_y = -\infty}^{\infty} \tilde{\bar{G}}^{ej}(k_{xm}, k_{ym}, z)\, e^{-j k_{xm} x}\, e^{-j k_{ym} y}\,}$$

with **Floquet wavenumbers**:

$$k_{xm} = k_{x0} - \frac{2\pi m_x}{d_x}, \qquad k_{ym} = k_{y0} - \frac{2\pi m_y}{d_y}$$

Each $(m_x, m_y)$ is a "Floquet mode" — a plane wave whose wavevector is the steering wavevector shifted by $2\pi/d$. Mode $(0,0)$ is the main beam; the rest are higher-order modes that may or may not propagate.

---

### 2.4 Method of Moments and Galerkin projection

You can't solve Maxwell's equations on the dipole metal exactly, so you make an **ansatz**: the current is a known shape (the basis function) times an unknown amplitude.

**Basis function** (entire-domain, sinusoidal):

$$b(x,y) = i(x) \cdot j_t(y)$$

$$i(x) = \frac{\sin\!\left(k_0\!\left(\tfrac{l}{2} - |x|\right)\right)}{\sin\!\left(k_0 \tfrac{l}{2}\right)} \quad \text{(textbook current on a thin dipole)}$$

$$j_t(y) = \frac{\mathrm{rect}_{[-w/2,\,w/2]}(y)}{w} \quad \text{(uniform across a narrow dipole)}$$

In the spectral domain its FT is $I(k_x) \cdot J_t(k_y)$ — the $I$ and $J_t$ that show up in the impedance formula:

$$B(k_x, k_y) = \frac{2 k_0 \left(\cos\!\frac{k_x l}{2} - \cos\!\frac{k_0 l}{2}\right)}{(k_0^2 - k_x^2)\sin\!\left(k_0 \tfrac{l}{2}\right)} \cdot \mathrm{sinc}\!\left(\frac{k_y w}{2}\right)$$

**Galerkin projection**: enforce the boundary condition ($\vec{e}_{\mathrm{tan}} = Z_{\mathrm{surf}} \cdot \vec{j}$ on the dipole, $0$ elsewhere) by *projecting* it onto the basis function itself (test function = basis function). That collapses an integral equation onto a single scalar equation:

$$Z_{\mathrm{MoM}} \cdot I_0 = (V_0 - I_0 Z_L) \cdot \langle \mathrm{rect}_{\mathrm{gap}}/\delta, \, b \rangle$$

For a small feed gap $\delta$ the right-hand inner product $\approx 1$, giving:

$$\boxed{\,I_0 = \frac{V_0}{Z_{\mathrm{MoM}} + Z_L}\,, \qquad Z_{\mathrm{in}} = Z_{\mathrm{MoM}}\,}$$

---

### 2.5 Active input impedance — the central object

For an isolated antenna $Z_{\mathrm{in}}$ is one number. For an *infinite array*, every element is excited — so the impedance one element sees depends on what every other element is doing, including the scan angle. That's "active":

$$\boxed{\,Z_{\mathrm{in,act}} = -\frac{1}{d_x d_y} \sum_{m_x} \sum_{m_y} G_{xx}^{ej}(k_{xm}, k_{ym}) \, [I(k_{xm})]^2 \, [J_t(k_{ym})]^2\,}$$

**Why "active":** it's the impedance when *all* elements radiate together, not the passive (single-element) one. $Z_{\mathrm{in,act}}$ depends on $(k_{x0}, k_{y0})$ — the scan direction.

**Decomposition into radiation + reactance:**

$$Z_{\mathrm{in,act}} = Z_{00} + Z_{\mathrm{ho}}$$

$$Z_{00} = -\frac{G_{xx}^{ej}(k_{x0}, k_{y0}) [I(k_{x0})]^2 [J_t(k_{y0})]^2}{d_x d_y}$$

$$Z_{\mathrm{ho}} = -\frac{1}{d_x d_y}\sum_{(m_x, m_y) \neq (0,0)} G_{xx}^{ej}(k_{xm}, k_{ym}) [I(k_{xm})]^2 [J_t(k_{ym})]^2$$

- $Z_{00}$: fundamental Floquet mode → real part = power radiated into main beam.
- $Z_{\mathrm{ho}}$: sum over evanescent modes → mostly reactive, stored near-field energy.

Crucially: when a higher-order mode goes from evanescent to propagating (a grating lobe is born), $Z_{\mathrm{ho}}$ jumps from purely reactive to also having a real part. That extra real part *steals* radiated power from the main beam — and on the impedance side it can pull the resistance violently, killing the match. That's **scan blindness**.

---

### 2.6 Active reflection coefficient

$$\boxed{\,\Gamma_{\mathrm{act}}(k_{x0}, k_{y0}, f) = \frac{Z_{\mathrm{in,act}} - Z_L}{Z_{\mathrm{in,act}} + Z_L}\,}$$

You demand $|\Gamma_{\mathrm{act}}| < -10$ dB across all in-band frequencies and all in-spec scan angles. That's your matching spec. $Z_L = 100\,\Omega$ in your design.

---

### 2.7 Grating lobes & the grating-lobe diagram

A Floquet mode $(m_x, m_y)$ propagates if

$$k_{zm}^2 = k_0^2 - k_{xm}^2 - k_{ym}^2 \geq 0$$

i.e. if $(k_{xm}, k_{ym})$ sits inside a circle of radius $k_0$.

**Grating-lobe diagram**: in the $(k_x, k_y)$ plane, draw a disc of radius $k_0$ at every $(2\pi m_x/d_x, \, 2\pi m_y/d_y)$. The $(0,0)$ disc — shifted by your scan $(k_{x0}, k_{y0})$ — is the **visible region**. A grating lobe appears the moment any other disc overlaps the visible region.

**Sub-half-wave rule:** if $d_x < \lambda/2$ and $d_y < \lambda/2$ *at the highest frequency*, the discs are far enough apart that they never overlap, *for any scan angle*. You're grating-lobe-free everywhere.

**Closed form for the angle a grating lobe first appears:**

$$\sin\theta_{\mathrm{GL}} = \frac{\lambda}{d} - \sin\theta_{\mathrm{scan}} \qquad \text{(1D, principal plane)}$$

So to keep $\theta_{\mathrm{GL}}$ imaginary up to scan angle $\theta_0$:

$$\boxed{\,d_x \leq \frac{\lambda}{1 + \sin\theta_0}\,}$$

You picked $d = 4.8$ mm, $\lambda_{\min} \approx 9.67$ mm, so $\lambda/d > 2$ across the band → no grating lobes at any scan, period.

---

### 2.8 Scan blindness vs. Wood–Rayleigh anomaly

These are related but distinct:

- **Wood–Rayleigh anomaly**: the *frequency/angle* at which a higher-order Floquet mode just transitions from evanescent to propagating ($k_{zm} = 0$). Always exists if the lattice allows it.
- **Scan blindness**: a sharp null in the *active element pattern* (radiation pattern × matching) that occurs near the Wood anomaly because (a) energy redistributes between modes, and (b) $Z_{\mathrm{in,act}}$ can develop a strong resonance that drives $|\Gamma| \to 1$.

**Sub-half-wave lattice $\Rightarrow$ no Wood anomaly in the visible region $\Rightarrow$ no scan blindness.** Your figure on slide 11 (smooth roll-off, no spikes) is the visual proof.

---

### 2.9 Windowing approximation & active element pattern

Infinite arrays radiate a Dirac delta in angle (zero beamwidth, infinite power) — so to compute a real radiation pattern you assume each element carries the *infinite-array* current, but truncate to $N_x \times N_y$ elements. The far field factorises:

$$E_{\mathrm{total}}(\theta, \phi) = \underbrace{\bigl[\text{isolated element pattern}\bigr]}_{\text{geometry only}} \cdot I_0(k_{x0}, k_{y0}) \cdot AF_x \cdot AF_y$$

with array factors

$$AF_x = \sum_{n_x = 0}^{N_x - 1} e^{j (k_x - k_{x0}) n_x d_x}, \qquad AF_y = \sum_{n_y = 0}^{N_y - 1} e^{j (k_y - k_{y0}) n_y d_y}$$

The **active element pattern (AEP)** is what you'd measure for one element when all neighbours are present and terminated:

$$\mathrm{AEP}(\theta, \phi) \propto j k_{z0} \, \bar{\bar{G}}^{ej}(k_{x0}, k_{y0}) \, I(k_{x0}) \, J_t(k_{y0}) \, I_0(k_{x0}, k_{y0}) \, \frac{e^{-j k_0 r}}{2 \pi r}$$

where the scan-dependent piece is

$$I_0(k_{x0}, k_{y0}) = \frac{V_0}{Z_L + Z_{\mathrm{in,act}}(k_{x0}, k_{y0})}$$

The AEP has a built-in $\cos\theta$ (projected aperture) and additional roll-off from $G$ and the basis-FT, which is why scan loss exists at all.

---

### 2.10 HPBW and array sizing

For a uniformly excited rectangular array of length $L = N \cdot d$:

$$\boxed{\,\mathrm{HPBW} \approx \frac{0.886\,\lambda}{L} \quad (\text{broadside, in radians})\,}$$

For $\mathrm{HPBW} = 2°$ at the lowest frequency (largest $\lambda$): $\mathrm{HPBW} = 0.0349$ rad $\Rightarrow L = 0.886 \lambda / 0.0349$. Convert to $N = L/d$.

- **SLL (sidelobe level)** for uniform amplitude: $\sim -13$ dB.
- **Tapered amplitude** (Taylor, Hamming, etc.): lower SLL but wider beam. Standard trade-off.

---

## 3. The "why did you do X" questions — answered

These are the ones you *will* get from the panel. Memorize them.

### Why a printed dipole, not a patch?

- **Bandwidth.** A dipole is naturally broader-band than a patch (a high-Q resonator). 12% fractional bandwidth is hard for a patch.
- **Profile.** Single-layer, no cavity — fits aircraft skin.
- **Polarization.** A linear dipole is the simplest building block; circular polarization can be added with a second crossed dipole.

### Why $h = \lambda/5$, not the textbook $\lambda/4$?

- $\lambda/4$ maximises broadside gain at one frequency — but it's a *peak*, so it falls off on both sides $\Rightarrow$ narrow band.
- A shorter $h$ shifts the peak above the band, so within the band you sit on the lower flank, where the response is **flatter**. Trades a bit of peak gain for bandwidth.
- For 12% bandwidth, flatness wins.

### Why $d_x = d_y = 4.8$ mm (sub-half-wave)?

- $\lambda_{\min}$ in band $\approx 9.67$ mm at $31$ GHz. $\lambda_{\min}/2 = 4.84$ mm.
- $4.8 < \lambda_{\min}/2 \Rightarrow$ no grating-lobe disc overlap at any frequency, any scan.
- **Bonus**: no Wood anomaly in the visible region $\Rightarrow$ no scan blindness.
- You could go smaller — but smaller cells mean more elements for the same aperture, more cost, more T/R modules.

### Why dipole length $4.5$ mm (not the half-wave $5.1$ mm at $f_c$)?

- A textbook half-wave dipole at the centre frequency would be $5.1$ mm — *longer* than the cell ($4.8$ mm). Physically impossible.
- So you push it to $\sim 4.5$ mm, which puts the dipole's natural resonance just above $f_c$. That keeps $\mathrm{Re}(Z_{\mathrm{in}})$ inside the matching window across the whole band rather than collapsing toward the edges.

### Why $Z_L = 100\,\Omega$, not 50?

- The active resistance of this geometry sits around $100\,\Omega$ at broadside in mid-band. Picking $Z_L$ to match the *natural* $\mathrm{Re}(Z_{\mathrm{in,act}})$ maximises bandwidth.
- $50\,\Omega$ would force you back through a matching network, narrowing the band and adding loss.
- $100\,\Omega$ is also more compatible with balanced feed structures (twin-line, integrated baluns) than $50\,\Omega$.

### Why is the E-plane the bottleneck and not the H-plane?

This is the question students get wrong most often. The clean version:

- **H-plane**: scanning moves $k_{y0}$. The basis FT is in $x$ — so $[I(k_x)]^2$ doesn't see the scan motion. Only effect: a $1/\cos\theta$ in the wave impedance, which makes $\mathrm{Re}(Z)$ *grow* with scan. Mismatch is gentle.
- **E-plane**: scanning moves $k_{x0}$ — *directly into* $[I(k_x)]^2$, which falls off as $\cos^2\theta$. On top of that, the GF carries another $\cos\theta$. $\mathrm{Re}(Z)$ drops below the matching window first.
- Net result: E-plane scan limit $\approx 28°$, H-plane $\approx 37°$. The system limit = the worst plane = $28°$.

### Why $55 \times 55$ elements?

- $0.886\,\lambda/L = \mathrm{HPBW}$. Pick the worst case (largest $\lambda$ $\Rightarrow$ lowest $f$), or the centre $f$ and accept slight overshoot at the band edge.
- $55 \times 55$ gives $1.97°$ at $f_c$ (just under $2°$) and $2.09°$ at the lowest frequency (4% over).
- $58 \times 58$ hits $2°$ even at the absolute worst case if you want a strict guarantee.

---

## 4. The curveball questions — anticipate these

### "Why didn't you use a dielectric superstrate / cavity / EBG ground plane?"

Each of those buys narrowband matching tricks; EBG buys a high-impedance surface but at the cost of scan blindness creeping in. The image-theorem PEC reflector is the **simplest** structure that gives broadband forward radiation. You picked it because it's robust and the physics is fully analytic — no hidden surface waves.

### "How well does your infinite-array approximation hold for finite $N \times N$?"

For $N \approx 55$, the central elements behave like infinite-array elements to within a fraction of a dB. Edge effects (cylindrical waves from edges, spherical from corners) decay as $1/\sqrt{r}$ and $1/r$ and mostly affect a 2–3-element strip near the rim. Mitigation: **dummy edge elements** or **amplitude tapering**.

### "Have you accounted for mutual coupling?"

Yes — *fully*. Mutual coupling is exactly what the higher-order Floquet sum $Z_{\mathrm{ho}}$ captures. The reactance you see in your impedance plot **is** mutual coupling. The infinite-array MoM doesn't ignore coupling, it includes it exactly under the periodic assumption.

### "Why 25 modes per axis ($51 \times 51 = 2601$ modes total)?"

Convergence study from Assignment 5: at 25 modes per axis the impedance is converged to within 1%. The high-order tail decays as $\sim 1/m^2$ (evanescent decay times the basis-FT roll-off), so doubling to 50 changes nothing measurable.

### "What about losses — conductor, dielectric?"

Not modelled — your code uses PEC for both the dipole and the reflector, no dielectric substrate. In practice, copper at $30$ GHz has surface resistance $R_s \approx 0.05\,\Omega/\square$, which adds a few percent to $\mathrm{Re}(Z)$. It would shift the matching slightly but not undo the design. A real design would include a thin ($\sim 0.1$ mm) low-loss dielectric (Rogers RT/duroid, $\varepsilon_r \approx 2.2$) for printability.

### "Why not circular polarization for satcom?"

Right answer: most LEO/GEO satcom uses circular polarization to be insensitive to Faraday rotation in the ionosphere and to aircraft attitude. A real design would use **two orthogonal dipoles** fed in $90°$ quadrature, or a sequentially rotated subarray. This project asked for one polarization — a CP version is a straightforward extension.

### "Why a Galerkin projection — what's the alternative?"

Galerkin (test = basis) gives a **symmetric** matrix and is variationally optimal — small errors in the basis function give second-order errors in $Z$. Alternative is **point matching** (test = $\delta$-function), which is simpler but converges more slowly.

### "What if the satellite is at a low elevation and you need to scan to $60°$?"

Your design fails — the spec is honest at $28°$. To reach $60°$ you'd need to *redesign*: smaller lattice (to keep grating lobes out at higher scan), shorter dipole (to keep $\mathrm{Re}(Z)$ in the matching band over a wider scan range), maybe a lossy ground-plane choke or a thicker reflector cavity. This is *exactly* the kind of follow-up engineering a thesis defence wants you to acknowledge.

### "How would adding a dielectric slab change things?"

Three effects:

1. $k_0$ effectively becomes $k_0\sqrt{\varepsilon_r}$ in the slab, so dipole length shrinks by $\sqrt{\varepsilon_r}$.
2. **Surface waves** can be excited if $h >$ some threshold — those couple energy out at endfire and cause scan blindness even with $d < \lambda/2$.
3. The GF gets new reflection/transmission terms (Fresnel at the slab interfaces).

This is the "FSS on dielectric" generalisation in §3.4 of your notes.

### "What is scan blindness mechanistically?"

Two-sentence version: it's when a higher-order Floquet mode crosses the visible boundary, becoming a propagating plane wave that the array couples into strongly. Energy goes into that lobe instead of the main beam, *and* $Z_{\mathrm{in,act}}$ develops a resonance that drives $|\Gamma| \to 1$, so the array also reflects nearly all input power back to the feed.

### "Wood anomaly vs Rayleigh anomaly?"

Same physical event, different naming traditions. Wood (1902) noticed sharp anomalies in optical gratings; Rayleigh (1907) explained them as the wavelength at which a diffraction order grazes the surface ($k_z \to 0$). In array language: the Wood–Rayleigh frequency is where a grating-lobe Floquet circle just touches the visible region.

### "Why does adding a ground plane *broaden* the bandwidth?"

Two reasons:

1. The factor $\bigl(1 - e^{-j 2 k_z h}\bigr)$ adds a frequency-dependent constructive interference that, for a sub-quarter-wave $h$, has a *flat* (not peaked) response across a wide band.
2. It exponentially suppresses higher-order Floquet modes — those carry the inductive reactance that, in free-space arrays, ramps up steeply with frequency. Killing that reactive tail flattens the $Z_{\mathrm{in}}$-vs.-$f$ curve dramatically.

---

## 5. Keywords cheat sheet

| Keyword | When you say it |
|---|---|
| **Floquet theorem / Floquet mode** | Any time you describe the periodic GF or the impedance sum |
| **Spectral domain / plane-wave decomposition** | Justifying why image theorem is trivial |
| **Image theorem** | The single piece of new physics you added |
| **Galerkin projection** | When asked how you solved the integral equation |
| **Entire-domain basis function** | Your one sinusoidal basis function |
| **Active input impedance** | Always — never just "input impedance" for an array |
| **Active reflection coefficient** | Same |
| **Active element pattern (AEP)** | Radiation pattern of one element in the lattice |
| **Visible region / grating-lobe diagram** | Geometric arguments about scan limits |
| **Wood–Rayleigh anomaly** | Textbook name for the scan-blindness onset |
| **Sub-half-wave lattice** | Justification for no grating lobes / blindness |
| **Mutual coupling** | What $Z_{\mathrm{ho}}$ captures, mention it once |
| **Windowing approximation** | How you got far-field patterns from infinite-array currents |
| **Half-power beamwidth (HPBW)** | Q3 |
| **Taper efficiency / sidelobe level (SLL)** | If asked about uniform vs. tapered amplitude |

---

## 6. Five questions to drill on right now

If you can answer these five aloud, smoothly, you're ready:

1. *"Walk me through how the image theorem modifies the Green's function and why $h = \lambda/4$ is the broadside optimum."*
2. *"Your active impedance is a sum over Floquet modes. What does the $m = 0$ term physically represent, and what do the $m \neq 0$ terms represent?"*
3. *"What is scan blindness, where does it come from in the equations, and why doesn't your design suffer from it?"*
4. *"You picked the E-plane scan limit as the system limit. Explain why the E-plane is more constrained than the H-plane in this geometry."*
5. *"How does HPBW scale with frequency, and why is the lowest frequency the worst case for your array sizing?"*

Practice saying each answer out loud once. The defence is half about *what you say* and half about *how confidently you say it*.

---

## 7. Quick numerical reference card

| Quantity | Value |
|---|---|
| Frequency band | $27.5 - 31$ GHz |
| Centre frequency $f_c$ | $\sim 29.25$ GHz |
| Fractional bandwidth | $\sim 12\%$ |
| Lattice period $d_x = d_y$ | $4.8$ mm |
| Reflector height $h$ | $1.85$ mm $\approx \lambda_c / 5.5$ |
| Dipole length $l$ | $4.51$ mm |
| Feed reference $Z_L$ | $100\,\Omega$ |
| Worst-case broadside $|\Gamma_{\mathrm{act}}|$ | $-17.5$ dB |
| Spec margin at broadside | $7.5$ dB |
| E-plane scan limit | $\pm 28°$ |
| H-plane scan limit | $\pm 37°$ |
| D-plane scan limit | $\pm 43°$ |
| Modes per axis (truncation) | $25$ ($51 \times 51$ total) |
| Array size (recommended) | $55 \times 55 = 3025$ elements |
| HPBW @ $f_c$ for $55 \times 55$ | $1.97°$ |
