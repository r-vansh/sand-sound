# Physics of vibrating plates and membranes that produce Chladni figures (for real-time simulation)

Reading notes on reliability:
- "Cited Findings" are taken from the linked source. Items marked *(search snippet)* were seen only in a search-result summary, not in the full page.
- Equations from the Tuan et al. papers were extracted by an automated page reader that had to reconstruct garbled superscripts; the structure is reliable, but check exponents against the PDF before publishing them.
- "Inferences" are my own derivations or arithmetic from the cited formulas. They are not sourced statements.
- "Gaps" include well-known textbook numbers that I could not verify against a source in this session. They are labelled as unverified.

Symbols: w = transverse displacement, h = thickness, rho = density, E = Young's modulus, nu = Poisson ratio, D = flexural rigidity, L = side of square plate, a = radius (circular) or half-side (where stated), k = wavenumber, f = frequency, omega = 2*pi*f.

---

## 1. Governing equations and usable mode shapes (square free plate driven at centre; circular plates)

### Takeaway
A metal Chladni plate obeys the fourth-order Kirchhoff-Love bending equation (not the membrane wave equation), and a square plate with free edges has no closed-form modes. Two usable approximations exist: (a) the exact-in-the-limit Ritz expansion in free-beam functions u_m(x)u_n(y) +/- u_n(x)u_m(y), and (b) the much cheaper cosine (Neumann-Helmholtz) modes cos(n*pi*x/L)cos(m*pi*y/L) used by Tuan et al., which reproduce measured point-driven patterns well at mid/high frequency but give wrong low-mode frequencies.

### Cited Findings

**Plate equation**
- Isotropic Kirchhoff-Love plate: `D * del^4 w = q(x,t) - 2*rho*h * w_tt` where Wikipedia's h is the half-thickness and `D = 2 h^3 E / (3 (1 - nu^2))`. In full-thickness notation this is `D = E h^3 / (12 (1 - nu^2))` and `D del^4 w + rho h w_tt = q`. — [Wikipedia: Vibration of plates](https://en.wikipedia.org/wiki/Vibration_of_plates); full-thickness form of D confirmed in [Tuan et al. 2022, Entropy (PMC8870825)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8870825/)
- Time-harmonic form: `(del^4 - k^4) psi = 0`, which factorises as `(del^2 + k^2)(del^2 - k^2) psi = 0`. The first factor is the Helmholtz (membrane-like, propagating) part; the second is the evanescent part that only matters near edges and the driver. — [Tuan et al. 2022, Entropy](https://pmc.ncbi.nlm.nih.gov/articles/PMC8870825/)
- Dispersion of bending waves: `f(k) = C k^2 / (2 pi)` with `C = sqrt(D / (rho h))` (the page prints C = D/(rho h); the square root is required dimensionally and reproduces their number). For aluminium E = 70 GPa, nu = 0.33, rho = 2700 kg/m^3, h = 1 mm: theoretical `C/2pi = 0.248` m^2/s; fitted to experiment `C/2pi ~ 0.226` (square plate) and `~ 0.220` (triangle). — [Tuan et al. 2022, Entropy](https://pmc.ncbi.nlm.nih.gov/articles/PMC8870825/)
- Tuan et al. 2015 (J. Acoust. Soc. Am. 137(4), doi 10.1121/1.4916704) used a 240 mm square aluminium plate, modelled the driven plate with the inhomogeneous Helmholtz equation, and found the measured frequency-wavenumber relation agreed with Kirchhoff-Love theory; substituting the resonant wavenumbers into the response function reconstructed the observed patterns. *(search snippet)* — [search summary of the 2015 JASA paper via preprints.org citing it](https://www.preprints.org/manuscript/202401.2096/v1)

**Free-edge boundary conditions (the hard part)**
- On each straight free edge (n = normal, s = tangential): bending moment `W_nn + nu W_ss = 0` and effective shear `W_nnn + (2 - nu) W_nss = 0`; at each corner the twisting moment vanishes, `W_ns = 0`. "For rectangular geometries these cannot be determined explicitly." — [R. Porter, Univ. of Bristol, "Eigenfrequencies and eigenmodes for a thin rectangular elastic plate with free edges" (2017)](https://people.maths.bris.ac.uk/~marp/abstracts/free-edge.pdf)
- Same conditions in Cartesian form for the square (-1,1)^2, with the eigenproblem `del^4 w = lambda w`. Kirchhoff (1850) solved only the circular plate; the free rectangular plate "resisted attack" (Rayleigh) until Ritz (1909). — [Gander & Wanner, "From Euler, Ritz, and Galerkin to Modern Computing", SIAM Review 54(4), 2012](https://www.unige.ch/~gander/Preprints/Ritz.pdf)

**Ritz mode shapes for the free square plate (accurate, still cheap)**
- Ritz's basis functions: `w_mn = u_m(x) u_n(y) + u_m(y) u_n(x)` and `w_mn = u_m(x) u_n(y) - u_m(y) u_n(x)`, where u_m are the modes of a free-free bar on x in [-1, 1]. — [Gander & Wanner 2012](https://www.unige.ch/~gander/Preprints/Ritz.pdf)
- Free-bar functions (normalised on [-1,1]), exactly as in the paper's code:
  - `u_0 = 1/sqrt(2)`; `u_1 = sqrt(3/2) x`
  - m even: `u_m = (cosh(k_m) cos(k_m x) + cos(k_m) cosh(k_m x)) / sqrt(cosh^2 k_m + cos^2 k_m)`, with `tan k_m + tanh k_m = 0`, root near `m*pi/2 - pi/4`
  - m odd: `u_m = (sinh(k_m) sin(k_m x) + sin(k_m) sinh(k_m x)) / sqrt(sinh^2 k_m - sin^2 k_m)`, with `tan k_m - tanh k_m = 0`, root near `(m - 1/2) * pi/2`
  — [Gander & Wanner 2012](https://www.unige.ch/~gander/Preprints/Ritz.pdf)
- Full-length free-bar roots (K = 2 k_m): `K = 4.7300, 7.8532, 10.9956, 14.1372, ...` — [Gander & Wanner 2012](https://www.unige.ch/~gander/Preprints/Ritz.pdf)
- Ritz's simplification for larger even m: drop cos k_m against cosh k_m, giving a cosine plus an edge-localised cosh correction; roots "differ only slightly from m*pi/2 - pi/4". — [Gander & Wanner 2012](https://www.unige.ch/~gander/Preprints/Ritz.pdf)
- Ritz used nu = 0.225 and a 6-coefficient expansion per symmetry class; his eigenvalues for the first (doubly antisymmetric) class were 12.47, 379.85, 1579.79, 2887.06, 5969.67, 14204.92 (exact matrix: 12.49, 379.34, 1559.28, 2899.93, 5961.32, 14235.30). His computed pitches agree remarkably well with Chladni's measured pitches. — [Gander & Wanner 2012](https://www.unige.ch/~gander/Preprints/Ritz.pdf)
- Porter's modern Ritz recipe: expand `W = sum alpha_mn w_m(x/a) w_n(y/b)`, split into four symmetry classes (SS, SA, AS, AA = symmetric/antisymmetric in x and y), and solve four real symmetric eigenproblems; matrix elements are given in closed form (no numerical integration). Convergence is slow; N = 12..48 terms per direction. Fortran code accompanies the report. — [Porter 2017](https://people.maths.bris.ac.uk/~marp/abstracts/free-edge.pdf)

**Cosine (Neumann-Helmholtz) approximation used for point-driven plates**
- Tuan et al. approximate free edges by Neumann conditions (d psi/dn = 0), "as a first approximation". For 0 <= x,y <= L: `psi_nm = (2/L) cos(n pi x / L) cos(m pi y / L)`, `k_nm = (pi/L) sqrt(n^2 + m^2)`, n,m = 0,1,2,... — [Tuan et al. 2022, Entropy](https://pmc.ncbi.nlm.nih.gov/articles/PMC8870825/); same modes in [Tuan et al. 2018, Sci. Rep. (PMC6052176)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6052176/)
- Driven response (inhomogeneous Helmholtz): `(del^2 + k^2) Psi = A delta(r - r_s)`, solution `Psi(r; k) = sum_n a_n psi_n(r)`, `a_n = A psi_n*(r_s) / (k^2 - k_n^2)`; Green's function `G(r, r_s; k) = sum_n psi_n*(r_s) psi_n(r) / (k^2 - k_n^2)`. — [Tuan et al. 2022, Entropy](https://pmc.ncbi.nlm.nih.gov/articles/PMC8870825/)
- The Green's-function Helmholtz approach works well for standard square and triangular plates but may not extend easily to more complex boundary conditions. *(search snippet)* — [Preprints.org 202401.2096](https://www.preprints.org/manuscript/202401.2096/v1)

**Circular plates**
- Axisymmetric solutions are `W(r) = C1 J0(lambda r) + C2 I0(lambda r)` with `lambda^4 = beta omega^2`, `beta = rho h / D` (full thickness). Clamped edge: frequency equation `J0(lambda a) I1(lambda a) + I0(lambda a) J1(lambda a) = 0`; mode shape `J0(lambda_n r) - [J0(lambda_n a)/I0(lambda_n a)] I0(lambda_n r)`; `omega_n = lambda_n^2 sqrt(D/(rho h))`. — [Wikipedia: Vibration of plates](https://en.wikipedia.org/wiki/Vibration_of_plates)
- Kirchhoff's circular-plate solution has figures consisting only of concentric circles and radii (diameters). — [Gander & Wanner 2012, quoting Ritz](https://www.unige.ch/~gander/Preprints/Ritz.pdf)

**Membranes**
- A plate or membrane is sometimes "approximately described by the wave equation" `del^2 w = (1/c^2) w_tt` in lab reports. — [Lindemann & Vella, UIUC Phys 406 report, 2015](https://courses.physics.illinois.edu/phys406/sp2017/Student_Projects/Spring15/T_Lindemann_A_Vella_Physics_406_Final_Report_Sp15.pdf)
- A widely copied simulation recipe uses fixed-edge membrane modes `sin(m pi x) sin(n pi y)` with `f_mn = k sqrt(m^2 + n^2)`; its own README flags that "the membrane model replaces the bending physics of elastic plates", that edges are fixed rather than free, and that its frequency scale is abstract. — [ratwolfzero/chladni_figures README](https://github.com/ratwolfzero/chladni_figures)

### Inferences
- **Membrane vs plate is the key fidelity choice.** Membrane: f proportional to k (f ~ sqrt(m^2+n^2)). Plate: f proportional to k^2 (f ~ m^2+n^2). The cosine *shapes* are the same family; only the frequency-to-wavenumber map and the edge condition differ. For a metal plate use the plate map. A real stretched membrane or speaker-driven film uses the membrane map with fixed (sin) edges.
- **Why the cosine approximation works at all:** away from edges the evanescent (del^2 - k^2) part decays within about one wavelength, so the interior field is Helmholtz-like. It fails for the lowest modes, where "interior" does not exist (see section 2 numbers).
- **Code-ready cosine recipe (approximate, Tuan-style):** with plate coordinates x,y in [0,L], `W(x,y) = sum_{n,m} cos(n pi x_s/L) cos(m pi y_s/L) cos(n pi x/L) cos(m pi y/L) / (k^2 - k_nm^2 + i*2*k*kappa)`, `k = sqrt(2 pi f / C)`, `C = sqrt(D/(rho h))`. The damping term is my addition in k-space; Tuan 2018 puts damping in frequency (section 3).
- **Code-ready Ritz recipe (accurate):** precompute k_m roots, tabulate u_m(x) on a 1-D grid, build modes as small sums of u_m(x)u_n(y) products per symmetry class (eigenvectors of Porter's matrix), store as textures. Only this route gives correct low-mode shapes (for example curved nodal lines near corners).

### Gaps
- I could not open Leissa's NASA SP-160 or Leissa 1973 (J. Sound Vib. 31, 257-293) to copy the standard free-square-plate table. Unverified from memory: for nu = 0.3, `omega L^2 sqrt(rho h / D)` is about 13.47, 19.60, 24.27, 34.80 (double), 61.09, ... Treat as unverified; Porter's nu = 0.225 values in section 2 are sourced.
- Free-edge circular plate eigenvalues were not found in a fetched source. Unverified from memory: `lambda^2 = omega a^2 sqrt(rho h / D)` about 5.25 for (2 diameters, 0 circles), 9.08 for (0,1), 12.2 for (3,0), 20.5 for (1,1), 21.5 for (4,0) at nu ~ 0.33; clamped plate lowest value 10.22. General free-plate mode shape is `[J_m(lambda r) + c I_m(lambda r)] cos(m theta)`.
- Circular membrane values not sourced here: modes `J_m(j_mn r / a) cos(m theta)` with Bessel zeros j_01 = 2.405, j_11 = 3.832, j_21 = 5.136, j_02 = 5.520 and f proportional to j_mn (unverified in this session, standard tables).
- The Tuan et al. 2015 JASA paper itself (paywalled) was not read; its content is represented through the open 2018 and 2022 papers by the same group.

---

## 2. Frequency versus mode number; real numbers for lab plates

### Takeaway
For thin plates frequency scales with wavenumber squared, so pattern "cell count" along a side grows like sqrt(f), not f: quadrupling the frequency doubles the number of nodal lines across the plate. Chladni's law f = C(m + 2n)^p (p ~ 2 for flat circular plates) is the circular-plate version of this; f proportional to (m^2 + n^2) is the rectangular version.

### Cited Findings

**Laws**
- Chladni's law: `f = C (m + 2n)^p`, m = number of nodal diameters, n = number of nodal circles; C and p depend on the plate. For flat circular plates p is roughly 2; for cymbals, handbells and church bells p ranges from 1.4 to 2.4 and can vary between mode families of one object. — [Wikipedia: Chladni's law](https://en.wikipedia.org/wiki/Chladni%27s_law)
- Rayleigh noted f is proportional to (m + 2n)^2 for large f; Rossing (Am. J. Phys. 50, 1982) proposed the modified form `f = C_n (m + 2n)^p`, fitted to flat and non-flat plates including cymbals. *(search snippet of a 1996 ASA abstract by Wilbur & Rossing)* — [auditory.org ASA abstract 5pMU8](https://auditory.org/asamtgs/asa96haw/5pMU/5pMU8.html)
- Simply supported rectangular plate (exact): `W_mn = sin(m pi x / a) sin(n pi y / b)`, `omega_mn = pi^2 (m^2/a^2 + n^2/b^2) sqrt(D / (rho h))` (full thickness h). — [Wikipedia: Vibration of plates](https://en.wikipedia.org/wiki/Vibration_of_plates)
- Free-edge square, cosine approximation: `k_nm = (pi/L) sqrt(n^2 + m^2)` and `f = (C/2pi) k^2`, so `f_nm = (C/2pi)(pi/L)^2 (n^2 + m^2)`. — [Tuan et al. 2022, Entropy](https://pmc.ncbi.nlm.nih.gov/articles/PMC8870825/)
- Orthotropic (rolled brass) plate: `K^2 = (pi/a)^2 [(1 - delta) n1^2 + (1 + delta) n2^2]` with small symmetry-breaking delta = 0, 0.014, 0.022 depending on cut angle; brass constants E_x = 107.7 GPa, E_y = 126.5 GPa. This lifts the (n,m)/(m,n) degeneracy and visibly changes patterns. — [Tuan et al. 2018, Sci. Rep.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6052176/)

**Accurate free-square-plate eigenvalues (Ritz, nu = 0.225)**
- Dimensionless `lambda = rho h omega^2 a^4 / D` with a = HALF side (plate is -a < x < a). First five per symmetry class, N = 48 values:

  | Mode | SS | SA = AS | AA |
  |---|---|---|---|
  | 1 | 0 (rigid) | 0 (rigid) | 12.4552 |
  | 2 | 25.9900 | 80.9341 | 320.969 |
  | 3 | 35.6495 | 235.500 | 375.420 |
  | 4 | 269.437 | 730.513 | 1526.65 |
  | 5 | 876.747 | 1104.26 | 2687.15 |

  — [Porter 2017](https://people.maths.bris.ac.uk/~marp/abstracts/free-edge.pdf)
- Ritz's own lowest value 12.47 for the same class (the "+"-shaped two-nodal-line mode family in his first symmetry type). — [Gander & Wanner 2012](https://www.unige.ch/~gander/Preprints/Ritz.pdf)

**Lab hardware and observed numbers**
- Tuan group: square aluminium plate, side 320 mm, thickness 1 mm, driven 200-5000 Hz at 0.1 Hz resolution, silica sand grain size 0.3 mm; also a 365 mm equilateral triangle. — [Tuan et al. 2022, Entropy](https://pmc.ncbi.nlm.nih.gov/articles/PMC8870825/)
- Tuan group (brass): square, side 280 mm, thickness 0.8 mm, oscillator at the centre, plate fixed at the centre by a screw. — [Tuan et al. 2018, Sci. Rep.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6052176/)
- UIUC: PASCO WA-9607 kit, 24 cm x 24 cm square aluminium plate and 24 cm diameter circular plate on a PASCO SF-9324 mechanical wave driver (a speaker voice coil with a post); resonances studied from 200 to 4000 Hz; specific square-plate resonances examined at 802 Hz and 970 Hz. — [Lindemann & Vella, UIUC 2015](https://courses.physics.illinois.edu/phys406/sp2017/Student_Projects/Spring15/T_Lindemann_A_Vella_Physics_406_Final_Report_Sp15.pdf)
- PASCO manual suggests sweeping roughly 100 Hz to 7 kHz, with aluminium plates about 0.035 inch (0.9 mm) thick. *(search snippet)* — [PASCO WA-9607 manual](https://cdn.pasco.com/product_document/Chladni-Plates-Kit-WA-9607-Manual.pdf)
- University of Toronto demo lists example pattern frequencies 280, 750, 1440, 2150 Hz for its plate. *(search snippet)* — [U of T Chladni plate demo](https://ugdemos.physics.utoronto.ca/db/demos/chladni-plate)
- An Oklahoma junior-lab report lists PASCO-plate resonances from 174 Hz to 2700 Hz. *(search snippet; the PDF's certificate had expired so it was not read)* — [OU junior lab report](https://www.nhn.ou.edu/~johnson/Education/Juniorlab/Chladni/2002-Chladni_S2002.pdf)
- Exact frequencies vary plate to plate: "the frequencies depend upon the thickness, exact size and strength of the metal"; same patterns appear, but not at quite the same frequencies. *(search snippet)* — [Science First EG-54 manual](https://sciencefirst.com/wp-content/uploads/2017/05/EG-54-chlandi-plates-manual.pdf)

### Inferences (my arithmetic from the cited formulas)
- **Aluminium constant:** D = 70e9 * (1e-3)^3 / (12 * (1 - 0.33^2)) = 6.55 N m; rho h = 2.7 kg/m^2; C = sqrt(D/(rho h)) = 1.557 m^2/s; C/2pi = 0.248 (matches the cited value). C scales linearly with thickness h and with sqrt(E/(rho(1-nu^2))).
- **Bending wavelength for 1 mm aluminium:** k = sqrt(2 pi f / C). At 250 Hz: k = 31.8 rad/m, wavelength 19.8 cm. At 1 kHz: k = 63.5, wavelength 9.9 cm (nodal-line spacing about 4.9 cm). At 4 kHz: k = 127, wavelength 4.9 cm. Rule: wavelength = 2 pi sqrt(C / (2 pi f)), proportional to 1/sqrt(f) and to sqrt(h).
- **Pattern complexity:** nodal lines across a side ~ k L / pi = 2L / wavelength. For L = 0.32 m, 1 mm Al: about 3 at 250 Hz, 6.5 at 1 kHz, 13 at 4 kHz. Number of nodal cells ~ (kL/pi)^2, proportional to f.
- **Cosine-model frequencies for the 320 mm, 1 mm Al plate:** f_nm = (C/2pi)(pi/L)^2 (n^2+m^2) = 23.9 Hz * (n^2+m^2) with the theoretical constant, 21.8 Hz * (n^2+m^2) with the fitted one. The 200-5000 Hz experimental range therefore spans n^2+m^2 from about 9 to about 230 (indices up to about 15).
- **Modal density is constant with frequency:** the count of (n,m) with n^2+m^2 < R is about pi R / 4, so there is roughly one mode per (4/pi)*21.8 = 28 Hz for that plate at all frequencies (about 110 Hz per centre-drivable mode once the even-even rule of section 3 is applied). This is why higher frequency does not mean sparser resonances.
- **Accurate low modes from Porter's table:** omega = sqrt(lambda) * C / a^2 = 4 sqrt(lambda) C / L^2. For the 320 mm, 1 mm Al plate f = 9.68 Hz * sqrt(lambda): AA-1 34 Hz, SS-2 49 Hz, SS-3 58 Hz, SA-2 87 Hz, SA-3 149 Hz, SS-4 159 Hz, AA-2 173 Hz, AA-3 188 Hz, SA-4 262 Hz, SS-5 287 Hz, SA-5 322 Hz, AA-4 378 Hz, AA-5 502 Hz. (Uses nu = 0.225 eigenvalues with nu = 0.33 aluminium stiffness, so expect errors of a few percent.)
- **The cosine model is wrong at the bottom:** it predicts (1,1) at 44-48 Hz and (2,0)/(0,2) at 87-96 Hz, while the true first three modes sit at 34, 49 and 58 Hz. Ratio of the true first three (nu = 0.225): 1 : 1.44 : 1.69, strongly Poisson-ratio dependent. Use the cosine model only above roughly the tenth mode, or accept it as a flagged approximation.
- Scaling a plate: all frequencies scale as h / L^2 and as sqrt(E / (rho (1 - nu^2))); patterns at corresponding modes are identical.

### Gaps
- No fetched source gave a complete table of (pattern, frequency) pairs for a specific commercial square plate. The UIUC report shows them only in a figure. The Science First manual points to Fletcher & Rossing, "The Physics of Musical Instruments", Table 3.2, which I could not open.
- Rossing's Am. J. Phys. 50 (1982) "Chladni's law for vibrating plates" was not read directly. The constant C for a specific plate was not found. Unverified from memory: Fletcher & Rossing give for the free circular plate's lowest mode f = 0.2413 c_L h / a^2 with c_L = sqrt(E / (rho (1 - nu^2))); a search for this formula did not confirm it.
- Chladni's own frequency ratios for the square plate and Ritz's comparison table are reproduced only as an image (Fig. 5.5) in Gander & Wanner; numbers were not extractable.

---

## 3. Between resonances: driven-damped modal superposition, selection rules, pattern morphing

### Takeaway
A point-driven plate is a sum over all modes, each weighted by (mode value at the driver) x (a complex Lorentzian in frequency); at a resonance one term (or one degenerate group) dominates, and between resonances several neighbours mix, so the nodal pattern deforms continuously instead of switching. Centre driving can only excite modes that are non-zero at the centre, which on a square leaves only the fully symmetric family, explaining why all observed centre-driven patterns have 4-fold symmetry.

### Cited Findings
- Modal response with damping (Tuan 2018, Methods): `Psi(x,y; omega) = sum_n [Q (m_d/m_p) / (1 + alpha Xi)] * [A omega^2 psi_n*(x',y') / (omega_n^2 - (omega - i gamma)^2)] * psi_n(x,y)`, with `Xi = sum_n A omega^2 |psi_n(x',y')|^2 / (omega_n^2 - (omega - i gamma)^2)`. Here (x',y') is the driving point, m_d/m_p the driver-to-plate mass ratio, alpha a dimensionless coupling factor. — [Tuan et al. 2018, Sci. Rep.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6052176/)
- Explicit square-plate coefficients: `C_{n1,n2} = [4 Q (m_d/m_p) / (1 + alpha Xi)] * (omega/omega_0)^2 cos(n1 pi x'/a) cos(n2 pi y'/a) / { [(1-delta) n1^2 + (1+delta) n2^2]^2 - [(omega - i gamma)/omega_0]^2 }`, with `Psi = sum C_{n1,n2} cos(n1 pi x/a) cos(n2 pi y/a)`. — [Tuan et al. 2018, Sci. Rep.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6052176/)
- Fitted parameters in that paper: coupling alpha = 0.1, mass ratio m_d/m_p = 0.1, damping `gamma = 0.16 omega_0`, omega_0 printed as "104 s^-1". — [Tuan et al. 2018, Sci. Rep.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6052176/)
- The governing driven equation there is an inhomogeneous biharmonic equation with a point source and a back-reaction term: `(del^4 - rho h omega^2 / D) Psi = (m_d omega^2 / D) [Q - alpha (m_p/m_d) Psi(x',y')] delta(x-x') delta(y-y')`. — [Tuan et al. 2018, Sci. Rep.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6052176/)
- Driver loading shifts resonances: without coupling, resonances are the poles of G (free-plate eigenvalues k_n). With coupling the response is `Psi_alpha = A G / (1 + alpha G(r_s, r_s; k))` and resonances satisfy `1 + alpha G(r_s,r_s;k) = 0`, giving redshifts that grow with alpha. The shape of |Psi_alpha(r_s,r_s;k)| is independent of drive amplitude A. — [Tuan et al. 2022, Entropy](https://pmc.ncbi.nlm.nih.gov/articles/PMC8870825/)
- In some of the new modes created by coupling, the driving point itself becomes a nodal point. — [Tuan et al. 2018, Sci. Rep.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6052176/)
- Experimental selection-rule evidence: centre-driven circular plate patterns "all exhibited circular symmetry"; centre-driven square plate patterns all "exhibited 4-fold rotational symmetry"; an off-centre-driven circular plate showed only mirror symmetry about the axis through driver and centre, and behaved "vastly" differently. — [Lindemann & Vella, UIUC 2015](https://courses.physics.illinois.edu/phys406/sp2017/Student_Projects/Spring15/T_Lindemann_A_Vella_Physics_406_Final_Report_Sp15.pdf)
- Sand patterns appear for only a fraction of the resonances a spectrum analyser detects. "Hidden" resonances lie very close to stronger ones; near 802 Hz and 970 Hz the real (in-phase) and imaginary (quadrature) displacement maps of the square plate look completely different, showing two modes vibrating out of phase with each other at one drive frequency. — [Lindemann & Vella, UIUC 2015](https://courses.physics.illinois.edu/phys406/sp2017/Student_Projects/Spring15/T_Lindemann_A_Vella_Physics_406_Final_Report_Sp15.pdf)
- Resonant patterns at higher frequency are superpositions of degenerate or nearly degenerate eigenmodes ("mode mixing"), attributed to degeneracy plus damping. *(search snippet)* — [Preprints.org 202401.2096](https://www.preprints.org/manuscript/202401.2096/v1)
- Resonances of real plates are sharp: PASCO advises changing frequency gradually because "the resonances are very sharp" and can be missed. *(search snippet)* — [PASCO WA-9406 manual](https://www.manualslib.com/manual/3028222/Pasco-Wa-9406.html?page=2)
- Modal loss factors of bare metal plates are "very low", low enough that mounting points measurably add damping; loss factors are usually measured by the half-power bandwidth method. *(search snippets)* — [Southampton eprint, "A study of joint damping in metal plates", Acta Astronautica 2009](https://eprints.soton.ac.uk/65878/1/2009_-_Acta_Astronautica_-_A_Study_of_Joint_Damping_in_Metal_Plates.pdf)
- A simple simulation convention: weight each mode by a real Lorentzian `1 / ((f - f_mn)^2 + gamma^2)`; small gamma gives one dominant mode, large gamma lets several contribute at once. Its author flags uniform gamma as an approximation. — [ratwolfzero/chladni_figures README](https://github.com/ratwolfzero/chladni_figures)

### Inferences
- **Centre-drive selection rule (cosine model).** The coefficient contains cos(n1 pi x'/a) cos(n2 pi y'/a). At the centre x' = y' = a/2 this is cos(n1 pi/2) cos(n2 pi/2), which is zero unless n1 and n2 are both even. So only even-even modes are driven, with sign (-1)^((n1+n2)/2). Three quarters of all modes are excluded.
- **Degenerate pairs lock into the "+" combination.** (n,m) and (m,n) have the same frequency and the same centre amplitude, so they are always driven equally and in phase: the excited shape is `cos(n pi x/L) cos(m pi y/L) + cos(m pi x/L) cos(n pi y/L)`. The "-" combination is zero at the centre and cannot be driven there. This is the physical origin of the popular "cos cos +/- cos cos" formula: for centre drive on an isotropic plate only "+" occurs; "-" patterns need off-centre drive, a bow at the edge, or anisotropy (the delta term) that splits the pair.
- **Exact-symmetry statement.** In Porter's classes, only SS modes (even in x and y) are non-zero at the centre, and among those only the ones also symmetric under x<->y. E.g. of the two low SS modes, the "diagonal cross" mode (u_2(x) - u_2(y) type) has a node at the centre and is not driven, while the "ring" mode (u_2(x) + u_2(y) type) is. AA-1, the lowest mode of the plate, is never excited by a centre driver.
- **Circular plate, centre drive:** only m = 0 (no nodal diameters) modes are non-zero at r = 0 since J_m(0) = 0 for m >= 1; patterns are concentric rings, matching the UIUC observation. Nodal diameters require off-centre drive or edge bowing/touching.
- **Lorentzian and Q from the cited denominator.** |omega_n^2 - (omega - i gamma)^2|^2 ~ (omega_n^2 - omega^2)^2 + 4 gamma^2 omega^2 for gamma << omega: a standard damped-oscillator line of full width 2 gamma (rad/s), Q = omega_n / (2 gamma). Phase runs from 0 below resonance through 90 degrees at resonance to 180 above.
- **Numeric Q estimate from Tuan 2018.** For their brass plate omega_0 = sqrt(B/(rho h)) (pi/a)^2 evaluates to roughly 1.1e2 s^-1 with handbook brass values (E ~ 115 GPa, rho ~ 8500, h = 0.8 mm, a = 0.28 m), so the printed "104 s^-1" is plausibly literal (not 10^4), and the fetched "(pi/a)" in the omega_0 formula must be (pi/a)^2 on dimensional grounds. Then gamma ~ 17 s^-1, linewidth 2 gamma / 2 pi ~ 5 Hz at every frequency, Q ~ 190 at 1 kHz and ~ 750 at 4 kHz. This is a mounted, sand-loaded, driver-coupled plate; a freely suspended bare plate would be sharper.
- **Why patterns morph smoothly.** Mode spacing (about 110 Hz between centre-drivable modes on the 320 mm plate) is much larger than the linewidth (about 5 Hz), so exactly on resonance one group dominates by a factor ~ spacing/linewidth ~ 20. Halfway between two resonances the two neighbours contribute comparably with opposite signs (one is above its resonance, one below), plus a smooth background from all distant modes. Sweeping f therefore moves nodal lines continuously: lines of mode A bend, reconnect at avoided crossings, and become lines of mode B. Amplitude drops by 1-2 orders of magnitude off resonance, so real sand stops moving there (pattern "freezes") even though the mathematical nodal set keeps evolving.
- **Nodal "lines" are strictly lines only without damping.** With damping the response is complex, W = W_re + i W_im. The time-domain motion is zero only where both vanish, which generically happens at points. Near a resonance one part dominates and |W| has deep narrow valleys that act as nodal lines. Simulate with |W|^2 = W_re^2 + W_im^2 and treat minima as the sand-collecting set. (Supported by the UIUC real/imaginary maps differing.)
- **Minimal faithful single-frequency algorithm:**
  1. `k2 = 2*pi*f / C` (C from material and thickness).
  2. For even n,m up to N: `k2nm = (pi/L)^2 (n^2+m^2)`, `omega_nm = C*k2nm`, `s = (-1)^((n+m)/2)`.
  3. Complex weight `a = s * omega^2 / (omega_nm^2 - omega^2 + 2i*gamma*omega)`.
  4. `W(x,y) = sum a * cos(n pi x/L) cos(m pi y/L)`; display |W|.
  5. Optional driver loading: divide by `(1 + alpha * Xi)` where Xi is the same sum evaluated at the driver. This only rescales overall amplitude and shifts peaks; it does not change the pattern at a given f.
  Truncate N so that omega_NN is several times omega; the series converges slowly at the driver point itself (log singularity of the 2-D Helmholtz Green's function), which is physical for a point force.
- The (n,m) = (0,0) term is rigid-body piston motion of the whole plate. Keep it only if the plate is modelled as free and pushed by a force; for a plate bolted to a shaker post, the more accurate picture is prescribed displacement at the centre.

### Gaps
- No fetched source tabulates measured Q factors or loss factors for Chladni plates. Unverified from memory: internal loss factors of aluminium and steel are of order 1e-4 to 1e-3 (Q of thousands), while mounted lab plates with sand typically show Q of order 100. The Q ~ 190-750 figure above is my inference from one fitted gamma.
- Whether damping is better modelled as constant gamma (constant linewidth in Hz, as in Tuan 2018) or constant loss factor (linewidth proportional to f) was not established. Radiation damping and air loading rise with frequency; not researched here.
- I found no source that plots nodal-line evolution during a continuous frequency sweep for a square plate; the morphing description is derived from the modal sum. A related paper title, "Visualization of Chladni Patterns at Low-Frequency Resonant and Non-Resonant Flexural Modes of Vibration", appeared in search results but was not read.
- The D'Alessio paper on forced vibrations and Chladni figures (SN Applied Sciences 2021, doi 10.1007/s42452-020-04062-6) and the Entropy 2024 paper "Exploration of Resonant Modes for Circular and Polygonal Chladni Plates" (mdpi.com/1099-4300/26/3/264) could not be fetched (login redirect and HTTP 403).

---

## 4. Several simultaneous frequencies and broadband music

### Takeaway
I found no experimental paper on multi-tone or music-driven Chladni plates; what follows is derived from linearity. Because the plate is linear, each frequency component produces its own complex response field and the time-averaged squared displacement is the sum of the squared magnitudes of the components, so exact nodal lines generally disappear and are replaced by minima of a weighted sum of single-frequency patterns, dominated by whichever spectral components fall on resonances.

### Cited Findings
- The driven response is linear in the source: `Psi = A * G(r, r_s; k)`; the response shape is independent of drive amplitude. — [Tuan et al. 2022, Entropy](https://pmc.ncbi.nlm.nih.gov/articles/PMC8870825/)
- Even at a single drive frequency, two modes can be simultaneously excited out of phase, producing different real and imaginary displacement maps (802 Hz and 970 Hz on a 24 cm square plate). — [Lindemann & Vella, UIUC 2015](https://courses.physics.illinois.edu/phys406/sp2017/Student_Projects/Spring15/T_Lindemann_A_Vella_Physics_406_Final_Report_Sp15.pdf)
- "Some of these modes might combine, so that the plate simultaneously undergoes different types of motion." *(search snippet)* — [UCSB lecture demonstrations, Chladni plates](https://web.physics.ucsb.edu/~lecturedemonstrations/Composer/Pages/44.45.html)
- Clear patterns require driving at resonant frequencies; Chladni sand patterns appear for only a fraction even of the true resonances. — [Lindemann & Vella, UIUC 2015](https://courses.physics.illinois.edu/phys406/sp2017/Student_Projects/Spring15/T_Lindemann_A_Vella_Physics_406_Final_Report_Sp15.pdf)
- Waller's two 1939-1940 Proc. Phys. Soc. papers on free square plates are titled Part I "normal vibrating modes" and Part II "compounded normal modes", i.e. combination patterns from simultaneously excited modes were catalogued as a distinct class. *(search snippet of a reference list)* — [MDPI Entropy 26(3) 264 reference list](https://mdpi.com/1099-4300/26/3/264)

### Inferences (derivation, not sourced)
- **Superposition.** For drive `F(t) = sum_j F_j cos(omega_j t + phi_j)`, the displacement is `w(x,y,t) = sum_j Re[ F_j W(x,y; omega_j) e^{i(omega_j t + phi_j)} ]` where W(.; omega) is the single-frequency complex field of section 3.
- **Time average.** For distinct frequencies the cross terms oscillate at the difference frequencies and average to zero over times long compared with 1/|f_i - f_j|: `<w^2>(x,y) = (1/2) sum_j |F_j|^2 |W(x,y; omega_j)|^2`. Mean-square velocity and acceleration use the extra weights omega_j^2 and omega_j^4: `<a^2> = (1/2) sum_j omega_j^4 |F_j|^2 |W_j|^2`. Sand bouncing depends on acceleration (compare with g), so high-frequency components count more per unit displacement.
- **Do nodal lines survive?** A point is motionless only if every component vanishes there, i.e. on the intersection of the nodal sets of all components. Two different patterns' nodal lines cross at isolated points, so with two or more comparable components there are no nodal lines, only nodal points, plus valleys where <a^2> is small. Sand then drifts to the minima of the summed field; pattern contrast is reduced and lines become blobs or broken segments.
- **Winner-takes-most.** Because each |W_j|^2 carries a resonance factor of order Q^2 (1e4 or more), a component sitting on a resonance outweighs an equally loud off-resonance component by a large factor. With a few tones the visible pattern is essentially that of the best-matched resonant component; clean patterns survive when one component dominates.
- **Special cases.** (a) Two tones within one linewidth of the same resonance drive the same mode: the pattern is unchanged and its amplitude beats at f1 - f2. (b) Two tones on two resonances of similar strength: the sand shows the union-minimum pattern (points where both patterns' nodal lines cross), the likely physical counterpart of "compounded" figures. (c) Harmonic complexes (a musical note): harmonics at f0, 2f0, 3f0... have wavenumbers scaling as sqrt(j), and they hit plate resonances essentially at random, so one or two harmonics usually dominate.
- **Music / broadband.** Treat a short-time spectrum S(f, t). Mean-square response is `<a^2>(x,y,t) = integral S(f,t) omega^4 |W(x,y; f)|^2 df`. Since W is a sum of narrow Lorentzians, this collapses to a sum over modes: `<a^2> ~ sum_n S(f_n) * (pi f_n / (2 Q_n)) * (omega_n^4 / (4 gamma^2 omega_n^2)) * psi_n(r_s)^2 * psi_n(x,y)^2` plus cross terms between overlapping modes. Practical form for code: one damped oscillator per mode driven by the audio signal (below), which handles all of this, including transients, automatically.
- **Time-domain modal synthesis (best fit for real-time audio):** for each retained mode n keep a complex or two-state oscillator: `q_n'' + 2 gamma q_n' + omega_n^2 q_n = psi_n(r_s) F(t) / M`. Integrate with an exact two-pole resonator (biquad) per mode at audio rate or block rate, then `w(x,y,t) = sum q_n(t) psi_n(x,y)`. Use the running mean of q_n q_m (or just q_n^2 if cross terms are dropped) to build <w^2>. Ring-up/ring-down time of each mode is 1/gamma (about 60 ms for gamma = 17 s^-1), so the plate has "memory": patterns lag the music by tens of milliseconds and sustained notes build much stronger patterns than percussive hits.
- **Cross terms matter for nearby modes.** If two modes lie within a linewidth of each other they are driven coherently by the same spectral component and their cross term does not average out; keep the full `|sum_n a_n psi_n|^2` per frequency component, and only add intensities across different frequencies.
- **Expectation for real music on a real plate:** mostly weak, smeared, slowly reorganising patterns; clear figures only during sustained tones that land within a few Hz of a centre-drivable resonance. This is physically correct and visually less dramatic than single-tone demos; an honest simulation should show that (a Q or "sustain" knob is physically meaningful because it corresponds to gamma).

### Gaps
- No peer-reviewed or lab source was found that measures Chladni patterns under two-tone, multi-tone or music excitation. Everything quantitative in this section is derived from linear superposition, which is safe for small amplitudes.
- Nonlinear effects at large amplitude (amplitude-dependent frequency shifts, mode coupling, subharmonics) are documented for plates and gongs by Rossing & Fletcher (1983, "Nonlinear vibrations in plates and gongs", which appeared in search results but was not read) and would break the simple intensity sum; not quantified here.
- How sand actually migrates in a field with no true nodes (rate versus local acceleration, the 1 g bouncing threshold, fine-powder inversion to antinodes) is outside this note's scope and is assumed to be covered by the companion research on grain dynamics.
- The content of Waller's "compounded normal modes" paper was not read; the link to simultaneous-mode patterns is inferred from its title.

---

## 5. Published patterns and frequency tables for square plates

### Takeaway
The classical catalogues are Chladni's own drawings (1787, 1802), Ritz's 1909 computed figures and pitch table, and Mary Waller's 1939-1961 classification; the modern quantitative references are Leissa's tabulated eigenvalues and the Tuan et al. point-driven experiments, which are the closest match to a speaker- or shaker-driven simulation. I could reach numeric tables only for Ritz/Porter eigenvalues; the rest are located but not extracted.

### Cited Findings
- Chladni published the technique in 1787 in "Entdeckungen über die Theorie des Klanges" (bowing a sand-strewn plate); a digitised copy is at ETH (doi 10.3931/e-rara-4235). Hooke observed nodal patterns on glass plates in 1680; Galileo noted such patterns in 1638. Chladni demonstrated to Napoleon in 1808, leading to the Paris Academy prize and Sophie Germain's work. — [Wikipedia: Ernst Chladni](https://en.wikipedia.org/wiki/Ernst_Chladni); e-rara reference via [ratwolfzero README](https://github.com/ratwolfzero/chladni_figures)
- Since the 20th century the bow is commonly replaced by a loudspeaker or electronic driver over or under the plate for precise frequency control. — [Wikipedia: Ernst Chladni](https://en.wikipedia.org/wiki/Ernst_Chladni)
- Ritz (Ann. Phys. 18, 737-807, 1909) computed the free square plate figures by his variational method and compared patterns and pitches with Chladni's measurements; full multi-term calculations for the first 9 notes, single-term for the rest; agreement is "remarkable". Gander & Wanner reproduce his figure comparison (Fig. 5.4), pitch table (Fig. 5.5) and give Maple/MATLAB code to regenerate everything. — [Gander & Wanner 2012](https://www.unige.ch/~gander/Preprints/Ritz.pdf)
- Before Ritz: Wheatstone (1833) approximated figures with sines and cosines (Ritz called this a trick applicable only in special cases); König (1864) suggested figures contain only straight lines, which Ritz showed to be wrong; Rayleigh also wrote "On the calculation of Chladni's figures for a square plate". — [Gander & Wanner 2012](https://www.unige.ch/~gander/Preprints/Ritz.pdf)
- Porter (2017) gives benchmark eigenvalues and nodal-line plots for the first 5 modes in each of the 4 symmetry classes of the free square plate (nu = 0.225), plus the 2:1 rectangle, with downloadable Fortran. — [Porter 2017](https://people.maths.bris.ac.uk/~marp/abstracts/free-edge.pdf)
- Leissa, "Vibration of Plates" (NASA SP-160, 1969) and "The free vibration of rectangular plates" (J. Sound Vib. 31, 257-293, 1973) are the standard compilations covering all 21 combinations of free/simply-supported/clamped edges. — cited in [Porter 2017](https://people.maths.bris.ac.uk/~marp/abstracts/free-edge.pdf)
- Mary Waller: "Vibrations of free square plates: Part I. Normal vibrating modes" (Proc. Phys. Soc. 1939) and "Part II. Compounded normal modes" (1940). *(search snippet of a reference list)* — [MDPI Entropy 26(3) 264](https://mdpi.com/1099-4300/26/3/264). Her 1957 work gives a rule relating the number of symmetry classes of a plate's normal modes to its number of lines of symmetry. *(search snippet)* — [Jain & Samajdar, "Nodal portraits of quantum billiards", arXiv 1709.03650](https://arxiv.org/pdf/1709.03650)
- Tuan et al. series (point-driven, sand on plate, automated frequency scans):
  - 2015, J. Acoust. Soc. Am. 137(4), "Exploring the resonant vibration of thin plates": 240 mm square and equilateral-triangle aluminium plates, inhomogeneous Helmholtz model. *(search snippet)*
  - 2018, Sci. Rep. 8, 10844: orthotropic brass squares, 280 mm x 0.8 mm, symmetry breaking and driver coupling. — [PMC6052176](https://pmc.ncbi.nlm.nih.gov/articles/PMC6052176/)
  - 2022, Entropy: 320 mm x 1 mm aluminium square and 365 mm triangle, 200-5000 Hz; resonant modes coincide with maxima of information entropy of the response. — [PMC8870825](https://pmc.ncbi.nlm.nih.gov/articles/PMC8870825/)
  - 2024, Entropy 26(3) 264, "Exploration of Resonant Modes for Circular and Polygonal Chladni Plates". *(search snippet; page returned 403)* — [MDPI](https://mdpi.com/1099-4300/26/3/264)
- Student-level data with photographs: UIUC Phys 406 reports (24 cm PASCO plates; sand patterns overlaid on near-field acoustic holography displacement maps, which agree closely). — [Lindemann & Vella 2015](https://courses.physics.illinois.edu/phys406/sp2017/Student_Projects/Spring15/T_Lindemann_A_Vella_Physics_406_Final_Report_Sp15.pdf); also [A. Powers, UIUC 2014](https://courses.physics.illinois.edu/phys406/sp2017/Student_Projects/Spring14/Anna_Powers_Physics_406_Final_Report_Sp14.pdf) *(search snippet, not read)*
- UNSW hosts still images and movies of a centre-shaker-driven square aluminium plate (no frequencies listed on the page). — [UNSW Chladni patterns](https://www.animations.physics.unsw.edu.au/jw/chladni-patterns.html)

### Inferences
- For validating a simulation, the best targets are (a) Porter's nodal-line plots and eigenvalues for the lowest modes (tests the Ritz path), and (b) the Tuan 2018/2022 figures of experimental versus computed patterns at stated wavenumbers (tests the point-driven cosine path, including driver coupling).
- Chladni's and Waller's catalogues were produced by bowing an edge while touching chosen points, which selects modes by clamping nodes with the fingers. They therefore include many figures (anything with a nodal line through the centre, all "-" combinations, antisymmetric classes) that a centre-driven plate can never show. A centre-driven simulation should be compared with centre-driven experiments, not with Chladni's full plate of drawings.
- Many attractive "Chladni" images online come from the formula `cos(n pi x) cos(m pi y) - cos(m pi x) cos(n pi y) = 0`; by the selection rule in section 3 that "-" family is exactly the one a centre driver does not excite on an isotropic plate. For centre drive use "+", with even n and m (in corner-origin coordinates).

### Gaps
- Waller's book "Chladni Figures: A Study in Symmetry" (1961) and her Proc. Phys. Soc. papers were not accessible; her (m,n) classification diagrams and frequency-ratio tables are not reproduced here.
- Chladni's original frequency (pitch) table for the square plate was located only as an image in Gander & Wanner (Fig. 5.5) and in an 1800s Phil. Trans. article excerpt without values.
- Tuan et al. give resonant frequencies and wavenumbers only in plotted spectra, not numeric lists, in the open papers; the extraction found no per-pattern (n1,n2) index table.
- Leissa's tables were not opened (see section 1 gaps for unverified remembered values).
- Fletcher & Rossing, "The Physics of Musical Instruments" (Table 3.2 and pp. 73-80 on plates) was not accessible online; it is the reference lab manuals point to for theoretical modal frequencies.
