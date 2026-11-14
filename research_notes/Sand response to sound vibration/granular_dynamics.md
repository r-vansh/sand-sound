# Granular dynamics: how sand responds to a vibrating surface and to sound

Notation used throughout: f = drive frequency (Hz), omega = 2*pi*f, A = displacement amplitude of the surface, g = 9.81 m/s^2, Gamma = A*omega^2/g (dimensionless peak acceleration), D or d = grain diameter, H or h = layer depth, N = H/D = depth in grain diameters.

Reliability note: items under "Cited Findings" come from the linked source. Several primary papers (Melo 1994 PRL, Umbanhowar 1996 Nature, van Gerner 2010 PRE full text, Jaeger-Nagel-Behringer RMP 1996, Knight 1995 PRE, Dorrestijn 2007 PRL) were paywalled or only retrievable as unreadable PDFs, so they are represented through abstracts, arXiv versions, or reviews that cite them. Anything I am supplying from general physics knowledge or my own arithmetic is under "Inferences" and labelled as such.

## 1. When grains start to move, how they move on a Chladni plate, and how the pattern forms

### Takeaway
Grains lose contact with the surface and bounce only where the local peak acceleration exceeds gravity (Gamma = A*omega^2/g > 1); above that they hop randomly with hop size set by the local vibration amplitude, so they behave like random walkers whose diffusivity falls to zero at nodal lines, where they collect. Below Gamma = 1 grains do not bounce and (per van Gerner et al.) instead roll slowly toward the antinodes, so the classic pattern needs Gamma > 1 over most of the plate.

### Cited Findings
- Definition: Gamma = 4*pi^2*f^2*A/g. "For Gamma > 1, on each cycle the layer loses contact and later collides with the container." — [Melo, Umbanhowar & Swinney 1995, arXiv patt-sol/9507003](https://ar5iv.arxiv.org/html/patt-sol/9507003)
- For a vibrated bed, "the threshold for collective motions to appear corresponds to the acceleration of gravity", i.e. Gamma about 1 (Evesque & Rajchenbach, as cited). — [Pastor et al., arXiv cond-mat/0402040](https://ar5iv.labs.arxiv.org/html/cond-mat/0402040)
- Classic Chladni mechanism: a strongly vibrating horizontal plate makes sand grains bounce and, over time, diffuse across its surface; this produces the Chladni figure, with grains collecting along the nodal lines of the resonating plate. — [Devauchelle et al. 2026, "Thermodynamics of bouncing grains", arXiv 2606.05930](https://arxiv.org/abs/2606.05930)
- The bouncing-grain diffusion is "mostly comparable to classical molecular diffusion": a temperature can be defined for the bouncing grains, the system obeys the fluctuation-dissipation theorem, and at uniform "temperature" grains follow Maxwell-Boltzmann statistics. When the temperature (vibration strength) varies across the plate, Fick's law no longer applies and a different transport relation proposed for Chladni's experiment holds. — [Devauchelle et al. 2026, arXiv 2606.05930](https://arxiv.org/abs/2606.05930)
- Random-walk picture: grains are described as random walkers whose steps shorten as they near a nodal line and which gather where their diffusivity is weak (attributed to Grabec and later work; this is from a search-result summary of the papers, not a read of the full text). — [Devauchelle et al., arXiv 2606.05930](https://arxiv.org/abs/2606.05930); [Chladni simulation write-up citing the stochastic model](https://dev.to/ratwolf/fast-chladni-pattern-simulation-with-a-simplified-membrane-model-using-vectorization-and-stochastic-30gf)
- When the plate resonates, motion is large over most of the surface, causing sand to bounce and move about; only at or near the node is the sand stationary. — [UNSW, Chladni patterns for guitar plates](https://newt.phys.unsw.edu.au/music/guitar/guitarchladni.html)
- Faraday's 1831 account: sand or heavy particles "proceed to the lines of rest"; powders "indicate the quiescent parts of the plates" and "figure out what are called the nodal lines". — [Faraday 1831, Phil. Trans. 121 (OCR transcription)](https://jimclifford.ca/philosophical-transactions-ocr-1665-1869/articles/jstor-107936/)
- Formation time (simulation): on a 40 mm x 40 mm plate in its 2x2 mode, 80,000 heavy particles (d = 0.075 mm) collect at the nodal lines "after a few seconds". — [van Gerner et al., J. Fluid Mech. 2011, "Air-induced inverse Chladni patterns"](https://cambridge.org/core/journals/journal-of-fluid-mechanics/article/airinduced-inverse-chladni-patterns/2A1779B2D6EB029119AC8C640BAF4B5F)
- Below Gamma = 1 the behaviour inverts: when the plate acceleration does not exceed g, particles move (roll) to the antinodes irrespective of size, even without air. The authors say this explains why standard Chladni experiments (which drive hard) missed it. — [van Gerner et al., Phys. Rev. E 82, 012301 (2010), abstract](https://research.utwente.nl/en/publications/inversion-of-chladni-patterns-by-tuning-the-vibrational-accelerat/)
- Demo practice: patterns observable from about 100 Hz up, example resonances 280, 750, 1440, 2150 Hz on a standard demo plate. — [U. Toronto demo page](https://ugdemos.physics.utoronto.ca/db/demos/chladni-plate)
- Demo practice: sweep from about 10 Hz to a few hundred Hz, listen for the loudest sound (resonance), then tune in 1-5 Hz steps. — [Rutgers physics demos 3D40.60](https://demos.physics.rutgers.edu/index.php?d=676)
- Use only as much sand as needed to show the nodal lines; at some high-frequency modes too much sand makes the standing wave hard to excite. — [UCSB lecture demonstrations](https://web.physics.ucsb.edu/~lecturedemonstrations/Composer/Pages/44.45.html)
- Sand bounces off the plate edge in demos (a catch box is recommended). — [UMass lecture prep 3D40.30](https://physicslectureprep.umasscreate.net/?p=16930)

### Inferences
- Amplitude needed for lift-off (my arithmetic from Gamma = 1, A = g/(2*pi*f)^2): 20 Hz -> 0.62 mm; 50 Hz -> 99 um; 100 Hz -> 25 um; 440 Hz -> 1.3 um; 1 kHz -> 0.25 um; 5 kHz -> 10 nm. So at audio frequencies the plate motion that throws sand is far too small to see; the plate looks still while the sand dances.
- Simulation model implied by the sources: each grain at position x sees local amplitude A(x) = A0*|mode shape(x)|. If Gamma_local = A(x)*omega^2/g > 1 the grain takes a randomly directed hop whose length scales with A(x) (ballistic estimate, mine: launch speed about A*omega, hop height about (A*omega)^2/2g, horizontal hop of similar order times a small random tilt). If Gamma_local < 1 the grain stays put (or, per van Gerner, creeps slowly up-amplitude). Diffusivity D_eff rises steeply with A(x) and vanishes in a band around each node.
- Consequence for line width (inference, not measured in sources): the resting sand line is the band where Gamma_local < 1. Near a node A(x) grows about linearly with distance, so band half-width w is about (g / (A0*omega^2)) / |grad(mode shape)|, i.e. w is of order (wavelength/2*pi)/Gamma_max. Driving harder (larger Gamma_max) gives thinner, sharper lines; driving just above threshold gives broad, fuzzy bands. Higher modes (shorter plate wavelength) also give thinner lines.
- Formation is fast near antinodes and slow near nodes (hop size shrinks), so patterns appear in seconds as a rough figure and then sharpen more slowly: outer grains arrive quickly, stragglers near the line take longest. Because hops shrink as grains approach, grains arriving later pile onto earlier ones, giving slightly heaped ridges rather than a single-grain line (inference; consistent with Faraday's "heaps" language but no ridge-height measurement found).
- Too much sand loads the plate, shifts resonant frequencies and damps the mode (supported qualitatively by the UCSB note above).

### Gaps
- No source found giving measured sand-line widths or ridge heights on real Chladni plates.
- No source read in full giving the measured law D_eff(A) or hop-length distribution; the Devauchelle et al. paper likely has it but only its abstract was accessible.
- Formation time on real plates is only anecdotal ("a few seconds" is from a simulation of a 40 mm plate).
- Dorrestijn et al. (PRL 2007, micro-scale Chladni figures on cantilevers) could not be retrieved; not covered here.

## 2. Why fine powder goes to the antinodes (inverse Chladni patterns)

### Takeaway
Two distinct, experimentally supported mechanisms send particles to antinodes: (a) air (or liquid) flow driven by the vibrating plate drags light, fine particles to the antinodes (Faraday 1831; confirmed and refined numerically by van Gerner et al. 2011), and (b) at sub-g accelerations any particle, of any size, drifts to the antinodes even in vacuum because it never bounces (van Gerner et al. 2010). The crossover to air-dominated behaviour is roughly where the grain's terminal fall speed drops to the plate's velocity, about 10-20 um grains at 50 Hz and Gamma = 1.

### Cited Findings
- Chladni already saw that fine powder collects at the antinodes, where vibration amplitude is largest; Faraday showed this depends on air by evacuating the container. — [Aranson & Tsimring review, arXiv cond-mat/0507419](https://ar5iv.labs.arxiv.org/html/cond-mat/0507419)
- Faraday 1831: lycopodium "produced the effect very well"; "chemically pure and finely divided silica rivals lycopodium in the beauty of its arrangement at the vibrating parts of the plate, although the same substance in sand or heavy particles proceeds to the lines of rest". So material is irrelevant; particle size/weight decides. — [Faraday 1831 (OCR transcription)](https://jimclifford.ca/philosophical-transactions-ocr-1665-1869/articles/jstor-107936/)
- Faraday showed powders on a plate are affected by currents of air as well as directly by the plate's vibration, with heaps of powder at the points of greatest vibration. — [The Book of Science summary of Faraday 1831](https://sharpgiving.com/thebookofscience/items/p1831c.html)
- Faraday described acoustic streaming in 1831; Rayleigh gave the theoretical explanation in 1884. — [APS Physics 9, 46 (2016)](https://physics.aps.org/articles/v9/46)
- Air-induced mechanism in detail: very light particles are carried to the antinodes by air currents set up by the plate. The time-averaged Eulerian air velocity actually points toward the nodal lines; what moves the particles is the cycle-averaged Lagrangian velocity of a tracer moving with the air, which points toward the antinodes. — [van Gerner et al., JFM 2011](https://cambridge.org/core/journals/journal-of-fluid-mechanics/article/airinduced-inverse-chladni-patterns/2A1779B2D6EB029119AC8C640BAF4B5F)
- Numbers from that simulation: 40 mm x 40 mm plate, 2x2 mode, 80,000 particles of d = 0.075 mm. With density 20,000 kg/m^3 they form the normal pattern in a few seconds; with density 20 kg/m^3 they form the inverse pattern after about 4 s. — [van Gerner et al., JFM 2011](https://cambridge.org/core/journals/journal-of-fluid-mechanics/article/airinduced-inverse-chladni-patterns/2A1779B2D6EB029119AC8C640BAF4B5F)
- Air becomes important when the grain's free-fall terminal velocity (Stokes: v_t = rho*g*d^2/(18*mu)) becomes comparable to the plate velocity; this holds for "10-20 um particles on a plate vibrating with frequency 50 Hz" at acceleration amplitude g. — [Aranson & Tsimring review, arXiv cond-mat/0507419](https://ar5iv.labs.arxiv.org/html/cond-mat/0507419)
- Sub-g mechanism: when plate acceleration does not exceed g, particles go to the antinodes irrespective of size, even in the absence of air. — [van Gerner et al., PRE 82, 012301 (2010)](https://research.utwente.nl/en/publications/inversion-of-chladni-patterns-by-tuning-the-vibrational-accelerat/)
- In liquid the effect is strong for big particles too: 0.75 mm glass beads on a 5 cm silicon plate submerged in water at 9,575 Hz go to the antinodes; water drag stops them lifting off, and the cycle-averaged horizontal force (surface motion plus gravity as the plate flexes) points toward the antinodes. — [APS Physics 12, 54 (2019) on Latifi, Wijaya & Zhou, PRL 122, 184301](https://physics.aps.org/articles/v12/54)
- Micro-scale version: polystyrene beads in water above a ~1 mm polysilicon membrane collect at antinodes via acoustic streaming; slightly off resonance they circulate in coordinated loops ("farandole"). — [APS Physics 9, 46 (2016) on Vuillermet et al., PRL 116, 184501](https://physics.aps.org/articles/v9/46)

### Inferences
- Rule of thumb for a simulation (my synthesis): compute Stokes terminal velocity v_t and plate peak velocity v_p = A*omega = Gamma*g/omega. If v_t >> v_p the grain is "sand" (goes to nodes when Gamma > 1). If v_t <~ v_p it is "dust" (advected to antinodes). Check (my arithmetic): quartz, d = 15 um -> v_t about 0.017 m/s; plate at 50 Hz, Gamma = 1 -> v_p about 0.031 m/s, consistent with the review's 10-20 um figure. Because v_p falls as 1/f at fixed Gamma, the dust/sand crossover size shrinks at higher frequency (about d proportional to f^-1/2 at fixed Gamma).
- Typical sand (100-500 um quartz, v_t of order 1 m/s or more) is firmly in the node-seeking regime in air at any audio frequency; flour, lycopodium (~30 um, low density), talc and fine silica are antinode-seeking. A mixture separates: coarse grains draw the lines, dust forms small round clouds/heaps at antinode centres.
- Three regimes by drive level for ordinary sand: Gamma_max < 1 everywhere -> nothing visible, or very slow creep to antinodes; Gamma_max slightly > 1 -> slow, broad pattern; Gamma_max >> 1 -> fast, sharp nodal lines.

### Gaps
- Exact particle sizes, frequencies and drift speeds in the van Gerner 2010 PRE paper were not accessible (abstract only). The speed of sub-g rolling toward antinodes is unknown here; expect it to be much slower than bouncing transport (inference).
- No quantitative Stokes-number threshold was retrieved from the JFM paper itself.

## 3. Deeper vibrated layers: heaping, convection, standing-wave patterns, oscillons

### Takeaway
A layer several grains deep behaves as a single inelastic body that leaves the plate each cycle once Gamma > 1; with air present it heaps and convects (Faraday heaping) from about Gamma 1-1.2, and from about Gamma 2.4-2.6 it shows subharmonic (f/2) standing-wave patterns: squares at low frequency, stripes at high, hexagons near Gamma 4, flat-with-kinks near 4.5-5.7, f/4 patterns above 5.7, and disorder above about 7.6. Pattern wavelength is a few layer depths and follows lambda/H = 1.0 + 1.1*(f*sqrt(H/g))^-1.32.

### Cited Findings

Regime map (bronze spheres 0.15-0.18 mm, layer 1.2 mm = 7 grains deep, f = 67 Hz, 127 mm cell, evacuated to 0.1 Torr):
- Gamma < 1: layer rides with the plate. Gamma > 1: layer leaves the plate each cycle but stays flat until about Gamma 2.4. — [Melo, Umbanhowar & Swinney 1995](https://ar5iv.arxiv.org/html/patt-sol/9507003)
- Gamma about 2.4: parametric wave instability; standing waves oscillating at f/2 (each crest becomes a trough on the next drive cycle). Squares for f < 24 Hz, stripes for f > 40 Hz. Stripes pictured at Gamma = 3.3. — [Melo et al. 1995](https://ar5iv.arxiv.org/html/patt-sol/9507003)
- Gamma about 3.9-4.0: f/2 hexagons (the layer's flight time exceeds one drive period: period doubling of the centre-of-mass motion). — [Melo et al. 1995](https://ar5iv.arxiv.org/html/patt-sol/9507003)
- Gamma about 4.5: hexagons give way to flat domains oscillating at f/2 with opposite phase, separated by "kinks" (pictured at Gamma = 5.8). — [Aranson & Tsimring review](https://ar5iv.labs.arxiv.org/html/cond-mat/0507419); [Melo et al. 1995](https://ar5iv.arxiv.org/html/patt-sol/9507003)
- Gamma about 5.7: squares/stripes reappear oscillating at f/4 (pictured at 6.0); f/4 hexagons at about 7.4. — [Melo et al. 1995](https://ar5iv.arxiv.org/html/patt-sol/9507003)
- Gamma about 7.6 and above: spatially and temporally disordered; persists to Gamma = 14 (largest tested). — [Melo et al. 1995](https://ar5iv.arxiv.org/html/patt-sol/9507003)
- Robustness: the same sequence was seen for particle sizes 0.05-3 mm, densities 2.3-11.4 g/cm^3, restitution 0.5-0.95, depths 2-40 particles; patterns are essentially independent of container shape. — [Melo et al. 1995](https://ar5iv.arxiv.org/html/patt-sol/9507003); [Umbanhowar, Melo & Swinney, Physica A 249 (1998) abstract](https://ideas.repec.org/a/eee/phsmap/v249y1998i1p1-9.html)
- Full catalogue: stripes, squares, hexagons, spirals, triangles, targets, plus localized oscillons and kinks. — [Physica A 249 (1998) abstract](https://ideas.repec.org/a/eee/phsmap/v249y1998i1p1-9.html)

Onset and hysteresis:
- Wave onset at Gamma about 2.6 on increasing drive; waves persist down to Gamma = 2.2 on decreasing (hysteretic). Pattern range studied 2.2 < Gamma < 3.8. Waves vanish for layers thinner than about N = 2 grains. — [Umbanhowar & Swinney 2000, arXiv nlin/0001061](https://ar5iv.labs.arxiv.org/html/nlin/0001061)
- The flat layer "remains flat for Gamma < 2.4", nearly independent of frequency. — [Aranson & Tsimring review](https://ar5iv.labs.arxiv.org/html/cond-mat/0507419)

Wavelength and scaling (grain diameters 0.05-0.8 mm, depths N = 1.1 to 30.7, bronze/steel/titanium/lead):
- Dimensionless frequency f* = f*sqrt(H/g). Wavelength: lambda/H = 1.0 + 1.1*(f*)^(-1.32 +/- 0.03). — [Umbanhowar & Swinney 2000](https://ar5iv.labs.arxiv.org/html/nlin/0001061)
- Square-to-stripe transition at f_ss = 0.33*sqrt(g/H) (f*_ss about 0.33; theory by Bizon et al. gives 1/sqrt(8) = 0.35). — [Umbanhowar & Swinney 2000](https://ar5iv.labs.arxiv.org/html/nlin/0001061)
- Grain mobility transition: grains within the layer move freely relative to each other only when peak container velocity v = 2*pi*A*f exceeds v_gm about 3*sqrt(D*g); below that the layer moves as a near-rigid block. — [Umbanhowar & Swinney 2000](https://ar5iv.labs.arxiv.org/html/nlin/0001061)
- 1994 result: at a critical acceleration there is a well-defined transition from flat to standing waves at half the drive frequency; the dispersion relation resembles that of gravity waves in a fluid but with a cutoff tied to particle size. — [Melo, Umbanhowar & Swinney, PRL 72, 172 (1994)](https://link.aps.org/doi/10.1103/PhysRevLett.72.172)

Bouncing-layer timing (inelastic bouncing ball model):
- Flight time reaches one drive period T at Gamma = sqrt(1 + pi^2) (about 3.30); period doubling at Gamma = sqrt(4 + pi^2) (about 3.72); flight time reaches 2T at Gamma = sqrt(1 + 4*pi^2) (about 6.36). Experiment (0.5 mm glass beads, 110 Hz) matches the model up to Gamma about 2.7 and shows the bifurcation at Gamma = 4.8 +/- 0.1. — [Pastor et al., arXiv cond-mat/0402040](https://ar5iv.labs.arxiv.org/html/cond-mat/0402040)

Oscillons:
- Stable, two-dimensional localized excitations in a vibrated granular layer; they oscillate at f/2 (a peak one cycle, a crater the next) and assemble into "molecular" and "crystalline" structures. — [Umbanhowar, Melo & Swinney, Nature 382, 793 (1996)](https://www.nature.com/articles/382793a0)
- They occur near the square/stripe crossover frequency f* (about 45 Hz for the Melo 1995 conditions) inside the hysteretic region just below wave onset; typical layers are 10-30 particle diameters thick. — [Aranson & Tsimring review](https://ar5iv.labs.arxiv.org/html/cond-mat/0507419)
- Opposite-phase oscillons attract at short range and bind into pairs/chains; like-phase repel; an oscillon persists for hundreds of thousands of cycles. Original experiment: brass balls, in vacuum. — [Wikipedia, Oscillon](https://en.wikipedia.org/wiki/Oscillon)
- Oscillons reported for layers at least about 14 monolayers deep (secondary essay; not verified against the Nature paper). — [UCSD course essay on oscillons](https://guava.physics.ucsd.edu/~nigel/Courses/Web%20page%20569/Essays_Fall2008/files/Gempel.pdf)

Heaping and convection (with air present):
- Heaps form slightly above threshold: fine grains vibrated at 1.2 g are "slightly above the threshold value for heap formation" (from search-result summary of the paper). — [van Gerner et al., "Interplay of air and sand: Faraday heaping unravelled", PRE 76, 051305 (2007)](https://link.aps.org/doi/10.1103/PhysRevE.76.051305)
- Mechanism: heaps are "formed and stabilized by the airflow through the bed" while the gap between bed and vibrating bottom is growing (air is sucked in under the lifted bed and drags grains inward/up the slope), confirming the pressure-gradient mechanism of Thomas & Squires (1998). An initially flat surface is unstable to heaping, and many small heaps coarsen into one large heap. — [van Gerner et al. 2007, TU/e record](https://research.tue.nl/en/publications/the-interplay-of-air-and-sand-faraday-heaping-unravelled/)
- Heaping diminishes or disappears when interstitial air can escape (permeable walls, evacuation) or when particles are larger (Pak et al. 1995). — [Aranson & Tsimring review](https://ar5iv.labs.arxiv.org/html/cond-mat/0507419)
- Convection: grains circulate in rolls; walls alone can drive it via friction. With 0.5 mm glass beads at 110 Hz a toroidal convection roll exists below Gamma about 6; convective drift is about an order of magnitude slower than the peak per-cycle grain velocity; behaviour is qualitatively unchanged above 60 Hz; sand gives the same qualitative results; cohesion matters for beads much smaller than 0.5 mm. — [Pastor et al., arXiv cond-mat/0402040](https://ar5iv.labs.arxiv.org/html/cond-mat/0402040)
- The wave-pattern experiments above were done in vacuum specifically to remove air-driven heaping. — [Melo et al. 1995](https://ar5iv.arxiv.org/html/patt-sol/9507003)

### Inferences
- Worked wavelengths from lambda/H = 1 + 1.1*f*^-1.32 (my arithmetic): H = 1.2 mm at 67 Hz -> f* = 0.74, lambda about 3.2 mm (stripes); H = 1.2 mm at 30 Hz -> f* = 0.33, lambda about 6.9 mm (square/stripe boundary). For H = 1.2 mm, f_ss = 0.33*sqrt(g/H) about 30 Hz, consistent with the reported squares below 24 Hz and stripes above 40 Hz. For H = 5 mm, f_ss about 15 Hz. So patterns are millimetres to about a centimetre across and coarsen (longer wavelength) at lower frequency and in deeper layers; at high f* the wavelength saturates near the layer depth H.
- Grain mobility threshold for D = 0.17 mm: v_gm about 3*sqrt(D*g) about 0.12 m/s (my arithmetic).
- These wave patterns need low frequencies and large strokes: Gamma = 3 at 30 Hz needs A about 0.83 mm, at 67 Hz about 0.17 mm (my arithmetic). They are a shaker-table or subwoofer-cone phenomenon (roughly 10-120 Hz), not something a 1 kHz tone produces; at 1 kHz Gamma = 3 is only 0.75 um of stroke, and f* is so high that lambda collapses to about H, i.e. no visible pattern, just a fluidised shimmer.
- Visual signature of the f/2 response: under steady light the surface looks like a blurred, doubled lattice; under a strobe at f/2 it freezes into crisp ridges. For a simulation, flip the phase of the height field every drive period.
- In open air with ordinary sand, expect the air-driven behaviours (a single heap with surface avalanches down the slope and upward flow inside; convection rolls) to mix with or mask the clean vacuum patterns, especially for fine sand. Coarse sand (>~0.5 mm) in a shallow layer is closest to the vacuum results. This is a synthesis of the Pak/van Gerner air results and the Melo vacuum results, not a single measured statement.
- Simple regime switch for a simulation of a deep layer: Gamma < 1 static (slow compaction only); 1-2.4 flat but bouncing, heaping/convection if air and fine grains; 2.4-4 f/2 squares (f < f_ss) or stripes (f > f_ss); ~4-4.5 hexagons; ~4.5-5.7 flat with kinks; 5.7-7.6 f/4 squares/stripes/hexagons; > 7.6 disordered.

### Gaps
- Exact oscillon parameter window (Gamma, f, depth) from the 1996 Nature paper not retrieved; commonly quoted values of about Gamma 2.45 and about 25 Hz are from my memory and unverified here.
- No numerical heap angle, heap growth rate or coarsening time found for Faraday heaping.
- Convection onset threshold beyond "about Gamma = 1" and roll sizes/velocities in absolute units not found.
- The explicit 1994 dispersion relation was not accessible; the 2000 scaling law above supersedes it for practical use.

## 4. Static sand properties: repose angle, avalanches, compaction, fluidisation, slumping

### Takeaway
Dry sand holds slopes up to about 34 degrees, fails by surface avalanches when over-steepened, compacts very slowly (logarithmically) under gentle tapping from about 58% toward 64% packing, and under vibration its friction is effectively reduced so heaps and ridges relax toward flat; strong vibration (Gamma >~ 2) flattens a pile outright.

### Cited Findings
- Angle of repose for dry sand is about 34 degrees. — [Wikipedia, Angle of repose](https://www.wikipedia.org/wiki/Angle_of_repose)
- 35 degrees listed for dry irregular sand. — [Dartmouth, Practical Fluid Mechanics ch. 16 table](https://cushman.host.dartmouth.edu/books/PracticalFM/Chap16.pdf)
- Two angles: a critical (maximum stability) angle just before an avalanche and a lower rest angle just after; exceeding the critical angle triggers avalanching until a stable angle is re-established. Rounded smooth grains pile less steeply than rough interlocking ones; moisture adds cohesion and allows far steeper slopes. — [UIUC Physics Van, Sand Piles](https://van.physics.illinois.edu/qa/listing.php?id=680); [Wikipedia, Angle of repose](https://www.wikipedia.org/wiki/Angle_of_repose)
- Compaction by tapping: density rises monotonically with tap number and is best fit by an inverse-logarithm law rho(t) = rho_inf - d_rho_inf / [1 + B*ln(1 + t/tau)], with all four parameters depending on tap strength (Knight et al. 1995). — [Peng & Ohta, arXiv cond-mat/9707237](https://ar5iv.labs.arxiv.org/html/cond-mat/9707237); [Talbot et al., arXiv cond-mat/9910239](https://arxiv.org/pdf/cond-mat/9910239)
- Scale: volume fraction rises from about 0.58 toward about 0.64 (random close packing); more than 10^5 taps are needed to approach steady state. — [Ben-Naim, Knight, Nowak et al. review text](https://cnls-www.lanl.gov/People/ebn/pubs/sandrev/sandrev.tex); [Ben-Naim et al., arXiv cond-mat/9603150](https://ar5iv.labs.arxiv.org/html/cond-mat/9603150)
- Heap relaxation: a heap is metastable; under strong vibration it relaxes to a horizontally flat surface even when its slope is already below the angle of repose. Simple linear (Fickian) diffusion does not describe the relaxation; a nonlinear transport law is needed. — [Tsuji, Otsuki & Katsuragi, arXiv 1708.03752](https://ar5iv.labs.arxiv.org/html/1708.03752)
- Related DEM/lab work used strong vertical vibration, Gamma >= 2, to drive heap flow. — [Tsuji et al., arXiv 1902.02942](https://ar5iv.labs.arxiv.org/html/1902.02942)
- Jaeger, Liu & Nagel (1989): with added high-frequency vertical vibration the pile angle relaxes logarithmically in time over a wide interval; vibration reduces the friction of individual surface grains. (Reported via a theory paper surfaced in search; the URL attribution was not verified by reading the document.) — [Augsburg thesis/paper on vibrated avalanches](https://opus.bibliothek.uni-augsburg.de/opus4/files/42052/42052.pdf)
- A vibrated heap spreads with an early super-diffusive then a sub-diffusive regime. — [EPL 2025](https://epljournal.edpsciences.org/10.1209/0295-5075/adc151)
- Fluidisation threshold: grains in the layer move relative to each other only above peak plate velocity about 3*sqrt(D*g); below it the bed moves as a block. — [Umbanhowar & Swinney 2000](https://ar5iv.labs.arxiv.org/html/nlin/0001061)

### Inferences
- For a height-field simulation: enforce a maximum slope tan(theta_c) with theta_c about 34-36 degrees and relax to theta_r about 2-4 degrees lower when it fails (the "few degrees" gap between the two angles is from my background knowledge of the granular literature, not a retrieved number). Make the effective repose angle a decreasing function of local vibration strength: at Gamma_local = 0 full angle; as Gamma_local approaches 1 creep with ln(t)-like slow slumping; for Gamma_local > ~1-2 treat the sand as fluid-like (slope relaxes toward zero within seconds).
- This explains why Chladni ridges are low and rounded rather than sharp 34-degree prisms: the flanks of each ridge sit where the plate is still shaking slightly, so the flanks keep slumping; only the core at the node is truly at rest. When the drive stops, whatever shape exists is frozen (friction returns instantly).
- Gentle vibration below Gamma = 1 does not move sand sideways visibly but does settle and compact it (surface sinks by a few percent: 0.58 -> 0.64 packing is about a 9% volume reduction at most, my arithmetic), and smooths small roughness.
- Compacted sand that has been shaken is slightly harder to re-mobilise than freshly poured sand (inference from the packing results; no direct source).

### Gaps
- Jaeger, Nagel & Behringer RMP 1996 and Physics Today 1996 full texts were not accessible; the specific numbers they give for the two angles were not retrieved.
- No source found for heap relaxation at Gamma < 1 in quantitative form (rates, angle vs Gamma curve).

## 5. Can airborne sound alone move sand? Membranes, speakers, Kundt's tube, cymatics

### Takeaway
Airborne sound does not push sand grains directly at any normal listening level (direct acoustic forces only matter around 145 dB and for very light particles); what works is sound driving a light, compliant surface (drum skin, stretched membrane, thin plate, the speaker cone itself) into resonance so that the surface exceeds 1 g, after which the sand behaves exactly as on a Chladni plate. Fine light powders in a resonant tube (Kundt) are the one classic case where the air itself organises the particles.

### Cited Findings
- Jenny's voice experiment: quartz sand on a black drum membrane 60 cm in diameter, set vibrating by singing loudly through a cardboard pipe, gave symmetrical Chladni patterns; low tones gave simpler, clearer figures and higher tones more complex structures. (Secondary source with promotional leanings; the physical observation is plausible and consistent with membrane modes.) — [Cymatics source document](https://www.cymaticsource.com/_files/ugd/0e8d92_35728f46be7e420c97d4109ad8edc33d.pdf)
- Jenny's main method: powders, pastes and liquids on a metal plate connected to an oscillator with a wide frequency range; substances formed geometric structures characteristic of the frequency. He used lycopodium powder (fine, light, very responsive) as well as sand. — [Wikipedia, Cymatics](https://en.wikipedia.org/wiki/Cymatics); [UCSC, Resonant Waves (Leonardo 2020)](https://creativecoding.soe.ucsc.edu/pdfs/Grillotti_ResonantWaves_Leonardo2020.pdf)
- Physical reading of cymatics: the nodal pattern is set by the shape (and boundary conditions) of the vibrating body and which of its modes the drive frequency excites; only the arrangement of sand changes. The phenomena "are best explained by classical physics"; Jenny's "force field" interpretation and later mystical/therapeutic claims are not established physics (Wikipedia files the topic under pseudoscience). — [Wikipedia, Cymatics](https://en.wikipedia.org/wiki/Cymatics)
- Each mode has a specific frequency, so the spectrum of the signal determines which patterns appear. — [Wikipedia, Cymatics](https://en.wikipedia.org/wiki/Cymatics)
- Kundt's tube (1866): fine powder (cork dust, talc, lycopodium) in a horizontal tube with a standing sound wave collects at the displacement nodes, in piles spaced lambda/2 apart; the detailed powder motion is due to acoustic streaming in the boundary layer at the tube wall. Driven originally by a stroked metal rod, now by a loudspeaker with a sine generator; the effect appears at resonance, when the sound "suddenly gets much louder". — [Wikipedia, Kundt's tube](https://en.wikipedia.org/wiki/Kundt%27s_tube)
- Within the agitated zones the dust also forms fine transverse striations (ripples about 1 mm or less apart) attributed to small vortices; cork dust is "blown into small piles". — [Britannica, Sound: open tubes](https://www.britannica.com/science/sound-physics/Open-tubes); [U. Washington PHYS 536 lecture notes citing Carman 1955](https://courses.washington.edu/phys536/lectures/536-09-resonators.pdf)
- Direct acoustic forces: levitation of extremely light objects begins around 145 dB in a standing wave. — [Frontiers for Young Minds 2021](https://kids.frontiersin.org/articles/10.3389/frym.2021.678021/full)
- Acoustic levitation in air typically uses 30-100 kHz ultrasound. — [arXiv 1808.03862](https://arxiv.org/pdf/1808.03862v1)

### Inferences
- Why air alone fails for sand (my estimate): in a plane wave, air particle acceleration is a = omega*p/(rho*c). At 100 dB SPL (p about 2 Pa) and 100 Hz this is about 3 m/s^2 for the air itself, but a quartz grain is about 2000 times denser than air and couples only weakly, so its acceleration is orders of magnitude below g. A rigid tabletop barely moves either. Hence sand on a table next to a loud speaker does nothing.
- Why a membrane works: a thin skin has tiny mass per area, so the pressure p moves it directly: near (not at) resonance a = p/(surface density). For a 0.1 kg/m^2 film and 2 Pa (100 dB) that is about 20 m/s^2, i.e. Gamma about 2, before any resonant gain (my arithmetic, idealised). So a stretched film, balloon skin or drumhead over a speaker or bowl is the realistic "sound moves sand" setup, and it responds only near its own resonances.
- Sand directly on a speaker cone / on a tray glued to the cone: the cone moves as a rigid piston at low frequency. Amplitude for Gamma = 1 is 0.62 mm at 20 Hz, 0.1 mm at 50 Hz, 25 um at 100 Hz (from section 1), all within reach of a woofer. Expected appearance by frequency and level (synthesis of sections 1-4, not a single cited observation):
  - Low frequency (20-100 Hz), low volume (Gamma < 1): nothing, or slow settling/compaction and smoothing of the surface.
  - Low frequency, medium volume (Gamma 1-2.5): whole layer bounces and looks blurred; sand migrates down any tilt, drifts toward the centre or the rim depending on cone shape, heaps and convects; grains get thrown off the edge.
  - Low frequency, high volume (Gamma 2.5-5+, deep enough layer): millimetre-to-centimetre ripples, squares or stripes flickering at half the drive frequency; above that, chaotic splashing.
  - High frequency (several hundred Hz to kHz) on a membrane or plate, enough volume to pass Gamma = 1 at antinodes: Chladni-type figures with thin lines; pattern complexity rises with frequency (more nodal lines, closer together); each figure appears only in a narrow band around a resonance and the sand sits still between resonances.
  - High frequency, low volume: nothing moves, since Gamma < 1 everywhere (required amplitude is tiny but so is the membrane's off-resonance response).
- Music rather than a pure tone: broadband signals excite many modes at once, so nodal lines of different modes overlap and no location stays quiet; sand jitters everywhere and collects only weakly, pulsing with bass hits. Clean figures need a sustained tone at a resonance. (Inference from the mode/spectrum statement above; consistent with common demo experience but I found no formal study.)
- Kundt spacing example (my arithmetic, c = 343 m/s): 1 kHz -> piles every 17 cm; 3.4 kHz -> every 5 cm.

### Gaps
- No reliable source found stating the SPL required to move sand on a membrane or the SPL at which Kundt-tube dust moves; the numbers above are my estimates.
- No peer-reviewed description found of sand directly on a loudspeaker cone under music; what is written above for that case is inference.
- Jenny's "Cymatics" volumes themselves were not accessed; the details come from secondary sources, some promotional. His tonoscope is not described in the sources retrieved.

## 6. What it looks like in motion

### Takeaway
While driven, sand is in constant visible agitation wherever the surface exceeds 1 g (a shimmering, fizzing haze of hopping grains that is busiest at antinodes) and dead still at the nodes; a pattern is a property of the resonance, so it reorganises abruptly when the frequency moves to a new mode and freezes in place the instant the drive stops.

### Cited Findings
- At resonance the plate's motion becomes large over most of the surface and the sand bounces and moves about; only at or near a node is it stationary. — [UNSW, Chladni patterns](https://newt.phys.unsw.edu.au/music/guitar/guitarchladni.html)
- Grain jumps become more dispersed near local extrema (antinodes) and particles settle near nodal lines. — [Stochastic models for Chladni figures, Proc. Edinburgh Math. Soc.](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/B834E2C39426B19BA7D9257772D94CB1/S0013091515000139a.pdf/stochastic-models-for-chladni-figures.pdf)
- Resonances are sharp: demonstrators tune in 1-5 Hz steps to find a figure, and the plate gets audibly louder at resonance. — [Rutgers physics demos](https://demos.physics.rutgers.edu/index.php?d=676)
- Low tones give simple figures, high tones complex ones (membrane with sand, Jenny). — [Cymatics source document](https://www.cymaticsource.com/_files/ugd/0e8d92_35728f46be7e420c97d4109ad8edc33d.pdf)
- Slightly off resonance in fluid-coupled systems, particles circulate in coordinated closed loops/vortices rather than sitting still. — [APS Physics 9, 46](https://physics.aps.org/articles/v9/46); [APS Physics 12, 54](https://physics.aps.org/articles/v12/54)
- In deep layers the surface pattern alternates each drive cycle (f/2): peaks become craters and back; oscillons look like single isolated peaks/craters that persist indefinitely. — [Melo et al. 1995](https://ar5iv.arxiv.org/html/patt-sol/9507003); [Umbanhowar et al., Nature 1996](https://www.nature.com/articles/382793a0)
- Small heaps merge over time into fewer larger ones (coarsening). — [van Gerner et al. 2007](https://research.tue.nl/en/publications/the-interplay-of-air-and-sand-faraday-heaping-unravelled/)
- Fine dust under the same conditions gathers in round heaps at the antinodes and keeps circulating within them, driven by air currents. — [Faraday 1831 (OCR)](https://jimclifford.ca/philosophical-transactions-ocr-1665-1869/articles/jstor-107936/); [van Gerner et al., JFM 2011](https://cambridge.org/core/journals/journal-of-fluid-mechanics/article/airinduced-inverse-chladni-patterns/2A1779B2D6EB029119AC8C640BAF4B5F)

### Inferences
- Sequence when a tone is switched on at a resonance: (1) instant "boil" over all antinodal regions; (2) within roughly the first second the antinode centres clear, because hops are largest there; (3) over a few seconds grains stream outward and pile against the quiet bands, drawing rough lines; (4) lines keep thinning and sharpening for tens of seconds as stragglers random-walk in with ever smaller steps; (5) steady state: still lines with a faint fizz of individual grains on their flanks. Timescales scale inversely with drive level; step 4 is intrinsically slow because diffusivity vanishes at the node.
- Frequency sweep: between resonances amplitude collapses, Gamma drops below 1 and the old figure simply sits there, frozen. On reaching the next resonance the old lines now lie partly on new antinodes, so they are thrown apart at once and re-gather on the new lines. Visually this is a sudden scatter followed by re-condensation, not a smooth morph. Sand that has left the plate or reached the rim is not recovered, so figures get sparser over a long sweep.
- Sound off: with no hopping there is no transport; the figure persists unchanged (dry friction, slopes below repose). Nothing "relaxes back". Erasing requires re-spreading the sand or driving a different mode.
- Volume up at fixed frequency: lines get thinner and cleaner, grains near the lines jump higher, more grains are lost off the edges; at very high drive even the nodal lines are disturbed by plate tilt and nonlinearity and the sand is sprayed around (inference; qualitative demo experience, not sourced).
- Volume down slowly: bands widen slightly as stragglers stop earlier, and ridge flanks stop slumping, leaving a slightly broader, softer line.
- Deep layer on a low-frequency shaker reads visually as a liquid: flat shimmering surface, then ripples; turning the drive off freezes ripples only if caught mid-phase, and usually leaves a near-flat, slightly dimpled surface because the last cycles relax it.

### Gaps
- No source found with quantitative time-lapse measurements of Chladni figure formation or of the transition between two figures during a sweep; the sequence above is reasoned from the transport physics and standard demonstrations.
- No measured numbers for how many grains are lost off the plate or how high grains jump as a function of drive.
