# 3D on 4D, polytope boundaries: balls (diameter 1), cubes and the other Platonic solids on the boundary of a regular 4-polytope

**Dimension.** Under Dimensions this is **3D on 4D**: solid 3D pieces packed on the 3-dimensional
boundary of a 4D polytope, with no container around them. Problem ids and records use the space
name `boundary`.

P balls of diameter 1 (radius 1/2) or P unit cubes (edge 1) are packed on the boundary of a tesseract of
edge s. Since 2026-11-02 the other four Platonic solids of edge 1 (tetrahedron, octahedron, dodecahedron, icosahedron), freely rotated, can be packed there too, and all five on each of the other five 4-polytopes (section 4C: on the tesseract the cube stays the axis-aligned one of section 4.3, everywhere else it is a rotated cube). Since 2026-11-01 balls can also be packed on the boundary of the other five regular
4-polytopes, the 5-cell, 16-cell, 24-cell, 120-cell and 600-cell (section 4A). Also since 2026-11-01 balls and
axis-aligned cubes (of dimension d - 1) can be packed on the boundary of the d-cube, d = 5 to 10 (4D on 5D to 9D on 10D, section 4B). The score is the smallest s,
the polytope's edge. Since 2026-10-27 balls have diameter 1; before, radius 1, and every ball score was twice as large (cube scores are unchanged). The sections below are stated in the site's units; the certifier and the ball code may work in radius-1 units internally (a factor of 2, exact) and convert at their edges. Every record is certified by an exact checker
(`python/boundary/certify.py`, exact arithmetic only), so its validity is proven, not just
tested. Records: `data/boundary/records.json`, from the PolyPack research prototype (2026-10-24),
the PolyPack boundary search (2026-10-25) and exact constructions (section 4.6). They are the tesseract's. The five other containers, and the Platonic solids other than the aligned cube, have
problems and lower bounds but no records (sections 4A.7 and 4C.6); their test packings are certified by the
general exact certifiers (`python/boundary/general.py`, section 4A; `python/boundary/platonic.py`, section 4C).

Labels: [EST] established and sourced, [CALC] computed or proved here, [SPEC] speculation.

## 1. The problem in plain words

A tesseract is a 4D cube. Its boundary is made of 8 solid cubes (the cells), glued face to face,
the way the surface of an ordinary cube is made of 6 squares. This boundary is a closed curved
3D space. It is flat inside every cell and across every square face (2-face) where two cells
meet. All its curvature sits on the edges: only 3 cells meet around an edge, so the angle around
it is 3 x 90 = 270 degrees instead of 360 [EST: Miller, Pak 2008, sec. 1: every point off the
faces of dimension at most d - 2 of a d-dimensional polytope boundary is flat; CALC: three cells
with 90-degree dihedral angles meet at each edge].

A piece is a rigid ball or cube that lives in this space. It may cross a 2-face into the next
cell, lying flat across it like a ribbon laid over the edge of a box. No edge or vertex of the
tesseract may pass through its interior. Pieces may touch but not overlap. The task: find the
smallest tesseract that holds P pieces.

There are 32 problems with records: P = 1..16 for balls and for cubes. Since 2026-11-01 the other containers have problems too (sections 4A and 4C): 816 on the six regular 4-polytopes in all.

## 2. Conventions

These are the conventions of the engine's contract (`packages/boundary/src/types.ts`).

- **Tesseract.** The cube [-s/2, s/2]^4 in R^4. Write h = s/2.
- **Cells.** Cell `2 i + (sigma > 0)` is the cube x_i = sigma h, for i = 0..3 and sigma = -1 or +1.
  So cell 0 is x_0 = -h, cell 1 is x_0 = +h, cell 2 is x_1 = -h, and so on up to cell 7, x_3 = +h.
  Cells 2i and 2i + 1 are opposite; every other pair of cells shares a 2-face.
- **Local frames.** A point of cell (i, sigma) has local coordinates u = the other three
  coordinates of R^4, in increasing order of their index. The cell is the cube [-h, h]^3 with its
  centre at the origin. Example: in cell 5 (x_2 = +h), u = (x_0, x_1, x_3).
- **Configuration.** `{s, pieces: [{cell, u, q?}]}`. `cell` is the anchor cell, the cell that
  holds the piece's centre; `u` is the centre in the anchor's local coordinates at this s, and it
  must lie in the closed cell. `q` is a unit quaternion (w, x, y, z) rotating a cube in the
  anchor's frame; absent means axis-aligned. Balls have no `q`.
- **Units.** Balls have radius 1/2 (diameter 1); cubes have edge 1 (the cube [-1/2, 1/2]^3 about its centre).
  Lengths scale with s, piece sizes do not.
- **Ids.** `boundary/tesseract/ball/<P>` and `boundary/tesseract/cube/<P>`, P without leading zeros.
- **Score.** s, the tesseract's edge. Smaller is better. Stored at full precision; shown at 6
  decimals (score_cert rounded to the nearest, ties to even, from its exact binary64 value).

Example, one ball: `{"s": 1.4142135623745096, "pieces": [{"cell": 0, "u": [0, 0, 0]}]}`.

## 3. The fold rule and what valid means

**Development.** A piece is placed by developing it from its centre c. Each point q of the piece
is reached by the straight segment from c to q. Where the segment meets a 2-face it continues
into the neighbouring cell, unfolded flat: in R^4 this is the rotation about the 2-face that takes
one cell's hyperplane onto the other's. This is the exponential map at c [EST: Miller, Pak 2008].
It is defined as long as no segment from c meets an edge first.

**Validity** (touching allowed):

1. each centre lies in its closed anchor cell;
2. no point of the 1-skeleton (the edges and vertices) lies in the interior of a developed piece;
3. each development is injective on the piece's interior (a piece does not overlap itself);
4. the interiors of distinct pieces are disjoint.

Given (2), conditions (3) and (4) together say: cut every developed piece into its parts, one per
cell and path class it reaches; in every cell, all parts (two parts of one piece included) have
disjoint interiors.

**Strict validity**: the closed pieces avoid the 1-skeleton, are embedded, and are pairwise
disjoint. Every record is certified strictly valid at its stored s (section 4.4).

**Edges** [CALC, research note sec. 2.2]. No angle condition is needed at an edge. A convex piece
fills at most 180 degrees around a boundary point, and every edge has 270 degrees.

**Vertices** [SPEC]. Whether a piece touching a vertex of the tesseract needs a further condition
(on its tangent cone there) is argued in the research note, not proved. The engine treats a vertex
like any edge point: inside a piece is invalid, on its boundary is allowed. The certified records
are strictly valid, so no piece touches any edge or vertex. The open question does not affect
them.

**nofold.** The simpler rule "every piece lies inside one closed cell" is not a separate track
here. Its best value is given with every record as a reference (section 6).

## 4. Certification

The certifier reads every number as an exact rational (a binary64 float is a dyadic rational, so
nothing is rounded) and decides every inequality in exact integer arithmetic: all coordinates are
put over one common denominator. It shares no code with the TypeScript engine.

### 4.1 Developments as lattice words

Every 2-face of the tesseract is a right angle in R^4. So when a straight path crosses a 2-face,
the next cell unfolds as the adjacent cube. All the cells a path crosses therefore unfold into the
cubes of the lattice 2h Z^3, with the anchor cell at the origin [CALC; the research prototype's
`tess.py`, re-derived]. Every edge of an unfolded cube lies on a lattice line (a line where two
lattice coordinates are odd multiples of h), and it is an edge of the tesseract.

A straight developed path crosses lattice planes y_k = (2n + 1) h. The ordered list of planes
crossed is a **word**. Along a straight segment each coordinate is monotone, so a word crosses the
planes of one axis all in one direction. Each word has an end cell and an exact development map
x = A y + h b from the lattice into R^4 (A a signed permutation, b an integer vector), built by
composing the unfoldings. Both depend only on the word, not on where the path starts.

The rotation used at each crossing: in cell (i, sigma), lattice direction tau e_k is the R^4
direction t e_j given by A. The path leaves through the 2-face x_i = sigma h, x_j = t h and
enters cell (j, t). The unfolding rotates in the plane of e_i and e_j about that 2-face, taking
sigma e_i (the old normal) to t e_j (the new normal) and t e_j (straight on) to -sigma e_i
(inward in the new cell).

**Realised word.** Word w realises a path from the start u to a target point x2 when the image
y = A^T (x2 - h b) lies in w's closed lattice cube and the segment u -> y crosses w's planes at
parameters 0 <= t_1 <= t_2 <= ... <= 1. The path's length is |y - u|. Two rules make the test exact:

- every cell that contains x2 is tried (x2 may lie on a 2-face shared by two cells);
- a plane the segment is parallel to is not counted as crossed. The segment then never meets it,
  or runs inside it. A path inside a 2-face plane cannot cross another plane without passing
  through an edge, so it stays in one closed cell, and the word without that crossing realises it.

Every realised word is a real curve on the boundary from u to x2. Ties t_j = t_{j+1}, or an end
on a lattice line, are curves through an edge. They are kept: they only add candidates.

**Completeness** [CALC]. *Every straight path of length L < 2s is realised by a word with at most
2 crossings per lattice axis.* Proof: 3 crossings in one axis means crossing three consecutive
planes, for example y_k = h, 3h and 5h. The segment then spans at least 5h - h = 4h = 2s in that
coordinate, so L >= 2s. For s > 1/2 this covers every path of length <= 1, which is all the ball
test needs. The certifier enumerates all such words: 1957 per start cell. A test repeats two ball
records (P = 8 and 9) with up to 3 crossings per axis (41179 words per cell) and finds the same
paths shorter than 2s.

The word maps are checked exactly in the tests: each unfolded cube lands on its end cell, a word
and its one-step extension agree on the 2-face they share, and going round an edge gives the
expected quarter-turn holonomy.

### 4.2 Balls

[CALC] A packing of balls of diameter 1 is valid exactly when

- every centre is at distance >= 1/2 from the 1-skeleton;
- every two centres are >= 1 apart (intrinsic distance);
- every geodesic loop at a centre (a straight path from the centre back to itself that avoids the
  skeleton) is >= 1 long.

The proof is that of [SURFACE.md](SURFACE.md) section 2.3, one dimension up. It uses one fact
about shortest paths [EST: Miller, Pak 2008, Prop. 1.2: a shortest path between two points that
do not lie in a common facet has no warped point in its interior; for a d-dimensional boundary
the warped points are the faces of dimension at most d - 2, here the edges and vertices]. Points
that do lie in a common cell are joined by the straight segment in that cell, which is shortest
because it is the R^4 distance. Strict validity uses > for all three.

How each condition is checked:

- **Clearance**, closed form. For a centre u in its cell, the distance to the skeleton is
  min over pairs of axes a < b of sqrt((h - |u_a|)^2 + (h - |u_b|)^2). Proof: a shortest path to
  the skeleton is a straight developed segment ending on a lattice line. The odd multiple of h
  nearest to u_a (which lies in [-h, h]) is +-h, so no lattice line is closer than this minimum,
  and an edge of the anchor cube itself attains it. The certifier compares the squared value
  with 1/4.
- **Pairs.** Intrinsic distance is never below the R^4 distance. So only pairs at most 1 apart in
  R^4 are examined; for those, every realised word from one centre to the other must give
  squared length >= 1. By completeness and Prop. 1.2, a pair closer than 1 always shows up as a
  shorter realised word. A realised word through an edge is at least clearance + clearance >= 1
  long, so it never causes a false rejection.
- **Loops.** Every realised non-empty word from a centre back to itself must give length >= 1.
  [CALC] Every word that ends in its own start cell (264 per cell) crosses two axes twice. So any
  loop shorter than 2s is at least s sqrt 2 long, which is >= 1 whenever a ball fits at all. The
  check is kept anyway. No record has a loop shorter than 2s (`min_loop` is null).

[CALC] Clearance <= h sqrt 2 = s / sqrt 2 for every point, so a ball needs s >= 1/sqrt 2 > 1/2. The
word stage therefore always runs at s > 1/2, where completeness holds.

### 4.3 Axis-aligned cubes

The certifier handles axis-aligned cubes, which is all the records use: q absent or
+-(1, 0, 0, 0). A rotated cube is refused as out of scope, not judged.

The cube centred at u is the box u +- 1/2 in the lattice picture. Say it **crosses axis k** when
|u_k| + 1/2 > h, so that it pokes through the face u_k = +-h.

[CALC] *It is valid alone exactly when it crosses at most one axis.*

- Two crossed axes k and l: the anchor's own edge u_k = +-h, u_l = +-h runs through the open
  box (the point with the third coordinate equal to u_m is inside it). Invalid.
- At most one crossed axis: in each other axis the open box, and u, lie inside the open
  interval (-h, h). So every segment from u to a point of the open box crosses u_k = +-h inside
  the 2-face, never at its rim. The box lies in the anchor cube and the neighbouring cube; no
  lattice line meets its interior. If s < 1, every axis is crossed, so one crossed axis forces
  s >= 1, and then the slice beyond the face (depth at most 1/2) stays inside the neighbouring
  cube.

The parts are then exact boxes: the box clipped to the anchor cell, and the slice beyond the
crossed face mapped into the neighbour by the one-crossing word. The map is a signed permutation,
so the slice is a box there too. The two parts of one cube lie in different cells, so a cube
never overlaps itself. Two parts in one cell have disjoint interiors when some axis separates
them (one's upper end <= the other's lower end).

Strict mode counts a face touched from inside (|u_k| + 1/2 = h) as crossed and adds the flat
slice in the neighbour. Closed boxes must then be separated (<). This catches two cubes that
touch across a 2-face from different cells.

### 4.4 score_cert

Scaling the tesseract by lam > 1 about its centre, with the pieces kept at unit size, is the same
as shrinking every piece about its own centre by 1/lam. A shrunk convex piece lies, with its
boundary, in the interior of the original one. So every contact of a valid packing becomes a gap,
and a valid packing becomes strictly valid [CALC]. The rule:

  for k = 1, 2, ...: lam_k = 1 + k 1e-12; scale s and every coordinate of every centre by
  lam_k (binary64 products); score_cert = s lam_k for the first k at which these floats certify
  strictly valid.

score_float is the unscaled s. The record's configuration is the scaled one, at s = score_cert,
so what is stored is exactly what was certified. Every record certifies at k = 1, except balls
P = 14 (k = 3), whose research result was 2e-12 short of touching.

### 4.5 Checks around the certifier

- **Records mode.** `python3 python/boundary/certify.py data/boundary/records.json` checks each
  record's own fields (id against polytope, piece and P; status; score_cert equal to the
  configuration's s; display; 0 < score_float <= score_cert; the stored bound below score_cert
  and below the proven bound, exactly) and then certifies the configuration. Every failure,
  malformed input included, is reported, and the exit status is 1.
- **Engine check.** `packages/boundary/scripts/check-records.ts` runs the engine's float
  validator on every record at score_cert, at the touching tolerance 1e-9 and at tolerance 0.
  All 32 pass.
- **Cross-check.** `packages/boundary/scripts/gen-crosscheck.ts` writes the engine's verdicts,
  chain-based distances, loops, clearances and cube parts for 314 configurations: the engine's
  research fixtures, the records, the records shrunk by 1e-6, the adversarial cases and 200
  seeded random ones (balanced valid and invalid). The test compares them with the certifier.
  Verdicts agree on all 314 at tolerance 1e-9 and on 287 at tolerance 0. The other 27 have an
  exact margin within 1e-13 of a threshold (raw float fixtures and exact contacts), where a float
  verdict at tolerance 0 is decided by rounding. One such case: the research fixture for 9 balls
  at its own s is exactly valid, with margins of about 1e-16; the engine rejects it at tolerance 0
  by rounding and accepts it at 1e-9. Distances, loops, clearances and cube part boxes agree to
  1e-9.
- **Adversarial tests.** Balls overlapping across a 2-face; balls overlapping around an edge
  (the shortest path crosses two 2-faces); a ball too close to an edge; a cube containing an edge
  point; two cube parts overlapping in one cell; touching variants of each, valid but not
  strictly; malformed and out-of-scope input. A copy of every record shrunk by 1e-6 fails.

### 4.6 Exact constructions

Some layouts sit on exact values. Where the certifier proves the exact version valid at exactly
that s (touching allowed), the record stores it in the field `exact`, records mode re-checks it,
and the record's configuration is that construction rounded to binary64 and scaled as in 4.4.

- **Cubes, rational** [CALC]. A cube layout whose value is within 1e-5 of 5/4, 4/3 or 3/2 is
  scaled to that value, and every coordinate within 1e-6 of a multiple of 1/24 (which holds the
  1/8, 1/6 and 1/4 grids) is snapped onto it; the other coordinates have slack and are kept as
  binary64 values. The certifier then decides validity at exactly s = 5/4, 4/3 or 3/2 in
  rational arithmetic. This works for the search layouts of P = 11, 12 (5/4) and 13, 14 (4/3) and
  the research layouts of P = 9 (5/4) and 14, 16 (3/2); it fails for the research layouts of
  P = 10, 12, 13 and 15, whose contacts run through coordinates off the grid. Those problems
  take a better or equal exact value from another candidate.
- **Balls, in Q(sqrt 2, sqrt 3)** [CALC]. Two layouts have every centre coordinate equal to
  +-a with h - a = 1/(2 sqrt 2), so every clearance is exactly 1/2:
  P = 4, s = 1/2 + 1/(2 sqrt 2) = (2 + sqrt 2)/4, a = (2 - sqrt 2)/8 (the research layout), and
  P = 8, s = cos 15 deg = (sqrt 6 + sqrt 2)/4, a = sin 15 deg / 2 = (sqrt 6 - sqrt 2)/8 (the
  search layout, one ball per cell). The certifier runs the same tests (clearance, words, loops)
  on elements a + b sqrt 2 + c sqrt 3 + d sqrt 6 with rational a..d (`python/boundary/qfield.py`).
  Signs are exact: a binary64 evaluation decides only when it exceeds 1e-12 times the sum of the
  terms' magnitudes, far above its proven error (2e-15 of that sum); otherwise sign(p + q sqrt 3)
  with p, q in Q(sqrt 2) follows from the signs of p, q and p^2 - 3 q^2, and the same one level
  down. Both layouts are valid at exactly their s, with every clearance exactly 1/2 and the least
  centre distance exactly 1. Whether these s are optimal is open [SPEC]: no matching lower bound.
- **P + 1 minus a piece.** A subset of a valid packing is valid, so the record for P + 1 with its
  last piece removed is a candidate for P. Cubes P = 10 (from 11) and 15 (from 16) use it.

The builder (`python/boundary/build_records.py`) certifies every candidate (research, search,
exact versions, P + 1 minus a piece, and the research note's constructions) with the rule of
4.4 and keeps the one with the smallest score_cert (ties: exact, search, research, subset).

## 4A. Balls on the other regular 4-polytopes

Added 2026-11-01. The engine is general over convex 4-polytopes, so the same rule (sections 2 and 3) holds on the other five regular 4-polytopes.
Only the tesseract's lattice words (4.1) are specific to it; for the others the developments are found by chains of cells, and the
exact certifier is `python/boundary/general.py`. Only balls are covered here (the Platonic pieces are section 4C).

### 4A.1 The containers

Each container has unit edge and is centred at its vertex centroid, as the contract says (`packages/boundary/src/types.ts`); the score is
its edge s, and a ball has diameter 1. The vertex lists are generated in a fixed order (`packages/boundary/src/containers.ts`; the
exact certifier repeats them in exact arithmetic) from the classical coordinates:

| container | vertices | edges | cells | ridges | cell | cells at an edge | angle around an edge | exact field |
|---|---|---|---|---|---|---|---|---|
| 5-cell | 5 | 10 | 5 | 10 | tetrahedron | 3 | 211.586 | Q (simplex in R^5) |
| tesseract | 16 | 32 | 8 | 24 | cube | 3 | 270 | Q (lattice words, 4.1; or chains) |
| 16-cell | 8 | 24 | 16 | 32 | tetrahedron | 4 | 282.115 | Q |
| 24-cell | 24 | 96 | 24 | 96 | octahedron | 3 | 328.414 | Q |
| 120-cell | 600 | 1200 | 120 | 720 | dodecahedron | 3 | 349.695 | Q(sqrt 5) |
| 600-cell | 120 | 720 | 600 | 1200 | tetrahedron | 5 | 352.644 | Q(sqrt 5) |

The angle is in degrees: the cells' dihedral angles at an edge, added over the cells there [CALC; REPORT.md 2.5 of the research note]. All six
are short of a full turn and above 180 degrees, so no angle condition is needed at an edge (section 3). The nearer to 360, the flatter the
container: the 120-cell and 600-cell are nearly flat (deficits 10.3 and 7.4 degrees), the 5-cell and 16-cell sharp [SPEC: the nearly flat ones
should behave like nofold packings until P is large].

- **5-cell**: the simplex with vertices (1, 1, 1, -1/sqrt 5), (1, -1, -1, -1/sqrt 5), (-1, 1, -1, -1/sqrt 5), (-1, -1, 1, -1/sqrt 5),
  (0, 0, 0, 4/sqrt 5) divided by sqrt 8. 16-cell: +-e_i / sqrt 2. 24-cell: the permutations of (+-1, +-1, 0, 0) / sqrt 2.
  600-cell: (+-1, 0, 0, 0) and permutations, (+-1/2)^4 and the even permutations of (+-phi, +-1, +-1/phi, 0) / 2, times phi. 120-cell: Coxeter's
  coordinates (the permutations of (0, 0, +-2, +-2), (+-1, +-1, +-1, +-sqrt 5), (+-phi^-2, +-phi, +-phi, +-phi), (+-phi^-1, +-phi^-1, +-phi^-1,
  +-phi^2) and the even permutations of (0, +-phi^-2, +-1, +-phi^2), (0, +-phi^-1, +-phi, +-sqrt 5), (+-phi^-1, +-1, +-phi, +-2)) divided by 3 - sqrt 5.
- **Cells, frames, ridges**: as section 2 and `polytope.ts` say: cells in the order of their outward unit normals (`compareNormals`), a cell's origin
  at its centre, its frame the standard basis projected on its hyperplane and orthonormalised in order, the gluing across a ridge the rotation of
  R^4 in the plane of the two normals that fixes the ridge. `buildPolytope` builds all six (the 120-cell, with 600 vertices, in about 0.4 s).
- **Ids**: `boundary/<polytope>/ball/<P>` with polytope `5cell`, `16cell`, `24cell`, `120cell` or `600cell`, P = 1 to 24 (the tesseract keeps
  its 32 problems, P = 1 to 16, ball and cube). The Platonic solids, a cube included, are section 4C (on these five the id `cube` is a rotated cube). 152 problems with balls and the tesseract's cubes; 816 with the solids.

### 4A.2 Chains, windows and completeness

A ball is valid exactly when its centre is at distance at least 1/2 from the hinges, every two centres are at least 1 apart and every geodesic
loop at a centre is at least 1 long (4.2; the proof is dimension-free and uses only Prop. 1.2 of Miller and Pak). What changes is how the
shortest paths are found. On the tesseract every ridge is a right angle, so the unfolded cells form a lattice (4.1). In general the cells a
straight path crosses unfold, one by one, into the hyperplane of the anchor cell (each gluing is an isometry that fixes its ridge), and the path
becomes a straight segment of that development.

**Chains and windows** [CALC; `packages/boundary/src/chains.ts` is the float twin]. From a centre u in its anchor cell, a *chain* is a sequence of
ridges R_1, ..., R_k, each a ridge of the cell reached by the previous one, developed into the anchor's hyperplane. Its *gate* is G_k = R_k
(developed) clipped by the *window* W_{k-1}, and the window is W_k = the cone from u over G_k (W_0 = everything). Lemma: a ray from u in the
interior of W_k crosses R_1, ..., R_k in this order and runs inside each developed cell between consecutive crossings, because each developed
cell lies on the far side of the plane of the ridge it was entered through and a ray crosses a plane once. A point y of the last developed cell
is therefore joined to u by the segment of the development, a path of the boundary of length |y - u|, exactly when y lies in the closed W_k.

**Completeness** [CALC]. *Every straight path from u that avoids the hinges and has length at most the reach L is realised by a kept chain,*
where the enumerator keeps a chain exactly when
1. u lies strictly on the inner side of the plane of the ridge (seen from the developed cell it leaves);
2. the clipped gate has positive area;
3. the gate is at distance at most L from u.

Proof. Take such a path with endpoint y and let R_1, ..., R_k be the ridges it crosses (finitely many, at points p_1, ..., p_k of their relative
interiors, since the path avoids the hinges, which are the ridges' edges). (1): the path leaves the anchor through R_1 at p_1; u lies on the inner
side of the plane of R_1, not on it, unless u itself lies on R_1, in which case the path starts into the neighbouring cell and is found from that
cell as its anchor (4A.3: a centre on a ridge is developed from both cells); likewise at every later ridge. (2): p_j lies in the interior of
W_{j-1}: a ray in a face of W_{j-1} lies in a plane through u and an edge e of the gate G_{j-1}, so it meets the plane of G_{j-1} on the line
of e, hence on e, since G_{j-1} is convex and e is one of its edges; e lies on a hinge or on a face of W_{j-2}, and by induction the ray meets
a hinge, which it does not. So a neighbourhood of p_j in R_j lies in W_{j-1}, and G_j has positive area. (3): |p_j - u| <= length <= L, and the
distance from u to G_j is at most |p_j - u|. So every prefix of the chain is kept and the path's endpoint lies in the closed window of the
last node; y in that node's cell, as a point of the developed cell. Hence a pair of centres closer than 1 along the boundary shows as a kept
chain whose window holds the image of the other centre at distance below 1: take L = 1, the ball's diameter. The same for loops (the image
of u itself, at a positive distance) and for the hinges (a hinge point at distance below 1/2 lies in the closed window of a kept chain, as the
segment from u to it, which meets no hinge before it, is such a path). *
Ties are harmless: an image on the boundary of a window lies on a path through a hinge point, of length at least the sum of two
clearances, at least 1 once both centres are 1/2 from the hinges.

**The turn rule** [CALC]. A straight line that meets no edge sees an edge under exactly pi in the plane normal to it, so it can cross consecutive
ridges around an edge only while the cells passed between them have dihedral angles at that edge adding up to less than pi (at most one cube of
the tesseract, 90 degrees, two tetrahedra of the 5-cell, 70.5 each, etc.). The float enumerator cuts a branch by it. The exact enumerator does
not need the rule: past pi around an edge the window, exactly computed, is empty (a ray cannot be on the far side of the edge), so the same
chains are kept. The exact enumerator's node counts equal the float one's on 316 of the 317 balls of the cross-check (one differs by 1, a window the float code cut as
thinner than MIN_WINDOW_ANGLE).

**Finiteness** [SPEC]. That the enumeration stops for every input is not proved here. The open windows of the kept chains of one depth are
pairwise disjoint (a ray that avoids the hinges follows one chain) and every gate is within L, which keeps the count small in practice: at
most 399 chains per ball in the test packings (4A.7), at most 245 on the 5-cell, 16-cell, 24-cell and 600-cell. The enumeration is capped (200,000
nodes in the exact enumerator, 4096 chains in the float validator, 1024 on the server) and a capped enumeration is refused, never accepted, so
a failure of finiteness could cost time and never soundness.

### 4A.3 The exact certifier

`python/boundary/general.py` (called by `certify.py` for the five polytopes, and for the tesseract in its tests as a second, independent
method). No floating point decides anything: every number is a Fraction, or an element a + b sqrt 5 of Q(sqrt 5) with rational a and b, with its
sign decided exactly (`Q5`: integer arithmetic, one gcd per operation; the rational containers use b = 0). It shares no code with the
TypeScript engine and repeats only its conventions (cell order, frames: the numbers that say what a configuration's (cell, u) means).

**Field choices.**
- The 16-cell and 24-cell (and the tesseract): coordinates in Z (the tesseract's in (1/2)Z), squared edge 2 (the 16-cell, 24-cell) or 1; the squared edge is rational.
- The 5-cell: the simplex with vertices e_i - (1/5)(1, 1, 1, 1, 1) in R^5, squared edge 2. The 5-cell has no rational embedding in R^4 (its
  coordinates there need sqrt 5 and sqrt 2), but its geometry is intrinsic: gluing rotations, windows and distances are affine and metric
  statements about vectors of a 4-dimensional hyperplane and need only the Gram matrices, so a rational copy in R^5 serves. A configuration's
  (cell, u) is read through four vertices of the cell (below), which is an affine map and needs no embedding.
- The 120-cell and 600-cell: Q(sqrt 5)^4 with the coordinates of 4A.1 (not divided by the edge: the base polytope has squared edge
  14 - 6 sqrt 5 and (3 - sqrt 5)/2).
- *Base units.* All lengths are those of the base polytope of squared edge e0^2; a ball of diameter 1 on the polytope of edge s has diameter
  e0/s there, so every test compares a squared base length with D2 = e0^2 / s^2, an exact number of the field (s is a binary64, a dyadic rational).
  Nothing needs a square root. Reported margins are converted by s^2 / e0^2 and rounded down.

**Reading a configuration.** The engine's (cell, u) means the ambient point s * origin + sum_k u_k frame_k (floats; the cell's frame is part of the
engine's conventions). The certifier writes that point as an affine combination of four independent vertices of the cell with dyadic weights
(rounded to 2^-60, summing to exactly 1), and the exact centre is the same combination of the exact vertices: a point of the cell's hyperplane,
within 3.4e-16 s of the engine's on the 5-, 16- and 24-cell and within 2.9e-15 s on the 120- and 600-cell ([CALC] over the review's corpus of engine readings). The strict test therefore asks for room: a pair distance or a loop must exceed 1 by more than 2e-13 s, and a clearance must exceed 1/2 by more than 1e-13 s (`READING_SLACK` in `general.py`; the reach is widened to match), so that "strictly valid" binds the configuration as the engine reads it. A centre that falls outside its closed cell by less than 1e-12 of an edge is moved onto the facet planes it
violates (as the engine's `clampToCell` does) before the exact test "the centre lies in its closed cell"; a larger violation is an error. The
cells that contain a centre are those whose hyperplane holds it, decided exactly (a centre on a ridge lies in two).

**Charts.** A cell's hyperplane has the affine chart x = O + lambda B (O a vertex, B three edge vectors of the cell, lambda in the field^3); the
metric is the exact Gram matrix G = B B^T. The gluing across a ridge between cells with outward normals N, N' (equal norms a, k = N.N'/a) is
the rotation
    Rot x = x + (k - 1) [ (N.x) N / a + (E.x) E / ((1 - k^2) a) ] + [ (N.x) E - (E.x) N ] / a,   E = N' - k N,
which takes N to N', fixes the orthogonal complement and has no square root (the tests check it is an exact isometry that fixes the ridge's
vertices, for every container). A development is a composition of such maps in charts (3 x 3 matrices of the field).

**Decisions** (all exact): the centre lies in its closed cell; for each ball and each cell containing it, the chains of 4A.2 are enumerated
(Sutherland-Hodgman clipping of the gate by the window's faces, cones from cross products of chart vectors, distances from squared lengths
with the Gram matrix); the clearance is the least squared distance from u to the part of each developed edge inside the window (the foot of the
perpendicular, clamped to the window's interval); for every ball j after i and every cell containing j, the image of j in each chain node
ending in that cell is tested against the closed window and its squared distance against D2; loops likewise for i with the trivial images (the
root, and the centre itself seen across its own ridge) removed. Strict validity asks for more than 1/2, more than 1 and more than 1; plain
validity for at least.

**Cost.** Container construction 0.01 to 0.5 s (the 120-cell, 600-cell: 0.4 s, gluings on demand); a ball needs 80 to 400 kept chains;
a test packing of 24 balls on the 24-cell takes about 4 s to certify, of 3 balls on the 120-cell about 7 s.

### 4A.4 Lower bounds

Both bounds of section 5 hold on every container [CALC]; each is stored as a binary64 rounded down by a relative 1e-15
(`packages/boundary/src/bounds.ts`; the exact check is `general.exact_bounds`, which compares with enclosures of pi, sqrt 2 and sqrt 5 to 45 digits).

- **Volume.** The developed balls have disjoint interiors, so P pi / 6 <= V s^3, V the volume of the boundary at unit edge, the cells' number times
  the cell's volume: 5-cell 5 sqrt 2 / 12, 16-cell 4 sqrt 2 / 3, 24-cell 8 sqrt 2, 120-cell 450 + 210 sqrt 5, 600-cell 50 sqrt 2 (the
  tesseract: 8).
- **One ball.** The centre needs distance at least 1/2 from every hinge, and every point of a regular cell lies within rho of an edge of that cell,
  rho the distance from the cell's centre to its edges (the cell's midradius: sqrt 2 / 4 for the tetrahedron, 1/2 for the octahedron,
  sqrt 2 / 2 for the cube, phi^2 / 2 for the dodecahedron). Proof: a point p of the cell lies on a segment from the centre O to a point q of a
  facet F of the cell; q lies within the inradius r_F of the regular polygon F of an edge e of F, since the farthest point of a regular polygon from its
  boundary is its centre; the distance to the convex set e is convex, so d(p, e) <= (1 - t) d(O, e) + t d(q, e) <= max(rho, r_F) = rho,
  because rho^2 = h^2 + r_F^2, h the distance from O to F's plane. So s >= 1 / (2 rho): sqrt 2 on the 5-cell, 16-cell and 600-cell, 1 on the
  24-cell, 1/phi^2 = 0.381966 on the 120-cell (and 1 / sqrt 2 on the tesseract, as in section 5).

Balls at the centres of two adjacent cells are 2 r s apart along the straight path across their common facet (r the cell's inradius, the path
through the facet's centre; any other path is longer), so P <= (number of cells) balls, one at the centre of each of P cells, are valid from
s0 = max(1 / (2 r), 1 / (2 rho)): sqrt 6 (5-cell, 16-cell, 600-cell), sqrt 6 / 2 (24-cell), 0.449028 (120-cell). The tests use these.
They are not optimal, and nothing about optimality is claimed [SPEC].

### 4A.5 Checks around the certifier

- **The tesseract again.** `general.py` also builds the tesseract (balls only). On the research records for balls (P = 4, 8 and 13 in the test; P = 3, 5, 8, 9 and 14
  by hand) it gives the clearance, the least distance and "valid" that the lattice-word certifier gives, to all digits: two
  independent methods (windows and charts against lattice words) agree.
- **Conventions.** `fixtures/containers.json` (written by `packages/boundary/scripts/gen-containers.ts`) holds the engine's vertex lists and, for
  many cells of each container, their vertex ids, origin and frame; the test compares them with the certifier's own reading to 1e-14.
- **Cross-check.** `fixtures/general_crosscheck.json` (`gen-general-crosscheck.ts`) holds 160 seeded random configurations of one to three balls (80
  valid and 80 invalid at tolerance 1e-9, 20 + 20 on each of the 5-cell, 16-cell and 24-cell, 10 + 10 on the tesseract, 6 + 6 on the
  600-cell and 4 + 4 on the 120-cell; the second ball placed at the end of a geodesic from the first, so that most pairs cross ridges) with the
  engine's verdicts, geodesic distances, loops and clearances. All verdicts agree (those within 2e-8 of a threshold are not compared: they depend on
  the engine's tolerance), and every distance, loop and clearance below the reach agrees to 1e-9.
- **Adversarial cases.** A ball centred on an edge of each container; one ball 1e-9 above and below its container's one-ball bound; two balls
  at the centres of adjacent cells at s0 (1 +- 1e-9) on every container (the valid side has its least distance within 3e-9 of 1); touching
  variants (the 24-cell at s = 1 is valid but not strictly valid; eight balls in the tesseract at s = 1 likewise); a centre at the middle of a ridge (it lies in two cells, and is developed from both; a second ball 0.9 away across the ridge overlaps it); a centre close to an edge of
  the 5-cell (cone angle 211.6 degrees) with a geodesic loop of length 0.0940692594805195, which the engine and the exact enumeration both find; a
  centre outside its cell; malformed input.
- **Completeness.** With 1.5 times the diameter as reach, instead of 1, the pairs, loops and clearances of length at most 1 are the same (exactly, on
  four test packings and on an overlapping pair).

### 4A.6 The problems, the server and the kernel

- **Problems.** 152 in all (664 more with the Platonic solids, 4C.1): the tesseract's 32 (P = 1 to 16, ball and cube), and balls on the 5-cell, 16-cell, 24-cell, 120-cell and 600-cell, P = 1 to 24 each
  (120). A smaller range for a container was not needed: the 24-ball problem is the most expensive to check (about 0.3 s on the 120-cell),
  and the range follows the research note's proposal (P = 1 to 24). The 120 new problems have no research record (status `no_record`): the first valid entry is kept, and an
  entry on a problem without a record must be valid at the strict tolerance 1e-11 (section 12).
- **Server.** The registry's boundary domain (`packages/boundary/src/domain.ts`, `packages/spaces/src/domain.ts`) gives each problem its lower bound
  (4A.4). The rule of section 12 is unchanged but for the cells (0 to C - 1, C = 5, 8, 16, 24, 120 or 600) and P (at most 24 on the five new containers).
  Checks of balls on the 120-cell from P = 5 and on the 24-cell from P = 22 (estimated 50 ms or more, `boundaryIsHeavy`) wait in the heavy-check
  queue, one at a time server-wide, and are charged their estimated time. The work is bounded: at most 24 pieces, 1024 chains per piece, volume
  and in-cell edge tests first. [CALC: the validator took 4 to 18 ms for 24 balls on the 5-, 16-, 24- and 600-cell and 91 ms on the 120-cell, loose
  layouts; a tight one costs up to three times that.]
- **Kernel** (section 14). `boundaryKernel` works on every container: the chain metric (`ChainMetric`) gives the energy, its exact gradient in the
  displacements and in the size (checked against central differences), the weights, the focus on one piece and the piece's clearance; poses are
  drawn by developing a uniformly random target through the chains from the piece's anchor at a growing reach (`kernel-pose.ts`; the shortest
  geodesic, tested by walking it back to the target); the seeds are balls at the centres of P cells (valid from s0) or, with more balls than
  cells, a sequential random layout; the bound is 4A.4's. The chain metric costs 0.12 (5-cell), 0.18 (16-cell), 0.25 (24-cell), 0.15 (600-cell) and 0.8 (120-cell) ms per ball for an
  evaluation with its gradient (about 100 chains each, 2 to 7 microseconds a chain), against 0.5 microseconds for the tesseract's words, so searches
  here are hundreds to a thousand times slower than on the tesseract. The board's rounds are sized to about two seconds of it (`boundaryRound`, apps/web).

### 4A.7 Test packings

No records: the new problems ship with lower bounds and test packings only (`python/boundary/fixtures/general_packings.json`, 25 packings; P = 1, 2, 3, 5,
8, 12 on the 5-cell, 1, 2, 4, 6, 10, 16 on the 16-cell, 1, 2, 3, 6, 12, 24 on the 24-cell, 1, 2, 3, 5 on the 600-cell, 1, 2, 3 on the 120-cell). Each is a
random sequential layout drawn in by the engine's `finish` until the balls touch (`gen-general-packings.ts`), then scaled as in 4.4
(lam = 1 + k 1e-12; k = 1 to 3) to the first k at which `general.py` certifies them strictly valid with every margin (clearance over 1/2, distances and loops over 1) at least 1e-12: 1.0e-12 to 2.0e-12 in all 25. They are not good packings (no search ran) and no optimality is claimed.

| container | P | s (score_cert) | lower bound | least centre distance | least clearance |
|---|---|---|---|---|---|
| 5cell | 1 | 1.469179 | 1.414214 (one ball) | more than 1 | 0.500000 |
| 5cell | 2 | 1.809986 | 1.414214 (one ball) | 1.101234 | 0.500000 |
| 5cell | 3 | 2.157183 | 1.414214 (one ball) | 1.000000 | 0.530271 |
| 5cell | 5 | 2.486605 | 1.643949 (volume) | 1.000000 | 0.572142 |
| 5cell | 8 | 3.084217 | 1.922774 (volume) | 1.000000 | 0.573311 |
| 5cell | 12 | 3.535233 | 2.201027 (volume) | 1.000000 | 0.580173 |
| 16cell | 1 | 1.549082 | 1.414214 (one ball) | more than 1 | 0.500000 |
| 16cell | 2 | 1.591073 | 1.414214 (one ball) | 1.535322 | 0.500000 |
| 16cell | 4 | 1.641365 | 1.414214 (one ball) | 1.031387 | 0.500000 |
| 16cell | 6 | 1.930761 | 1.414214 (one ball) | 1.000000 | 0.585401 |
| 16cell | 10 | 2.155838 | 1.414214 (one ball) | 1.000000 | 0.523211 |
| 16cell | 16 | 2.826915 | 1.643949 (volume) | 1.000000 | 0.551189 |
| 24cell | 1 | 1.034038 | 1.000000 (one ball) | more than 1 | 0.500000 |
| 24cell | 2 | 1.177151 | 1.000000 (one ball) | more than 1 | 0.500000 |
| 24cell | 3 | 1.166934 | 1.000000 (one ball) | 1.555277 | 0.500000 |
| 24cell | 6 | 1.181813 | 1.000000 (one ball) | 1.000000 | 0.506162 |
| 24cell | 12 | 1.265138 | 1.000000 (one ball) | 1.000000 | 0.512942 |
| 24cell | 24 | 1.417640 | 1.035623 (volume) | 1.000000 | 0.511808 |
| 600cell | 1 | 1.484763 | 1.414214 (one ball) | more than 1 | 0.500000 |
| 600cell | 2 | 1.676617 | 1.414214 (one ball) | more than 1 | 0.500000 |
| 600cell | 3 | 1.671674 | 1.414214 (one ball) | more than 1 | 0.500000 |
| 600cell | 5 | 1.654522 | 1.414214 (one ball) | more than 1 | 0.500000 |
| 120cell | 1 | 0.396759 | 0.381966 (one ball) | more than 1 | 0.500000 |
| 120cell | 2 | 0.434681 | 0.381966 (one ball) | more than 1 | 0.500000 |
| 120cell | 3 | 0.422123 | 0.381966 (one ball) | more than 1 | 0.500000 |

### 4A.8 Limitations

- **Platonic pieces** (the tetrahedron, cube, octahedron, dodecahedron and icosahedron, rotated, on the six containers) were not built when 4A was written; they are section 4C. Its plan for them, in order, and where each
  part is: (1) the float validator's parts for these shapes: `validate.ts` `validatePieces`, `parts.ts`, `pieces.ts`; (2) a penalty model for rotated pieces in the search and the kernel: `search/solids.ts`, `kernel-solids.ts`
  (dim 6, with a swap); (3) the exact certifier's polytope-part test: `platonic.py`; (4) the field: `qmulti.py` (degree 4 only for some pieces: 4C.3); (5) test packings: 4C.6.
- The vertex rule (section 3) is argued, not proved, on every container.
- The ball centre is read through float frames (4A.3): the certified configuration is within 3.4e-16 s (5-, 16-, 24-cell) or 2.9e-15 s (120-, 600-cell) of the one the engine sees; the strict test keeps a slack of 1e-13 s to 2e-13 s for it (4A.3).
- The 120-cell and 600-cell are nearly flat. Their problems are open, and their bounds are those of 4A.4; nothing is known about their optimal packings.

## 4B. The d-cube: (d-1)D on dD, d = 5 to 10

Added 2026-11-01. Under Dimensions these are **4D on 5D**, **5D on 6D**, ..., **9D on 10D**: balls of diameter 1 and axis-aligned unit (d-1)-cubes on the
boundary of the d-cube [-s/2, s/2]^d. It generalises the tesseract (section 2 to 4) with 4 replaced by d, and nothing else: the d-cube is the one regular
polytope whose every ridge is a right angle in every dimension, so the tesseract's lattice words carry over to Z^(d-1) [EST: research note 4.3; CALC: re-derived
below]. The boundary has 2d facets (the cells), each a (d-1)-cube of edge s; the curvature sits on the (d-3)-faces (the hinges), around each of which
three cells meet at right angles, 270 degrees, in every d. The score is the edge s. Ids are `boundary/<d>cube/ball/<P>` and `boundary/<d>cube/cube/<P>`
(containers `5cube` to `10cube`). Code: `packages/boundary/src/ncube/`; exact certifier: `python/boundary/ncube.py`. Nothing is shown in Tighter, and
the site's viewers do not draw these containers yet.

### 4B.1 Conventions

- **Cells.** Cell `2 i + (sigma > 0)` is the facet x_i = sigma h, h = s/2, i = 0..d-1, 2d cells. Its local coordinates u are the other d - 1 coordinates in
  increasing order, so the cell is [-h, h]^n with n = d - 1 (the tesseract's convention for n = 3).
- **Configuration.** `{s, pieces: [{cell, u}]}`, u of n numbers. There is no `q`: a cube is axis-aligned, and a rotation of a cube by one of its symmetries
  changes nothing (the tesseract's engine accepts such a `q`; here it is refused as malformed, so there is one form of every packing). A rotation that is not
  a symmetry of the cube is out of scope, as on the tesseract (4.3): the problems are about balls and aligned cubes.
- **Units.** Balls have diameter 1, cubes edge 1; lengths scale with s except piece sizes.
- **Problems and P.** Both pieces on each container, P = 1 to 160 (5cube), 384 (6cube) and 400 (7cube to 10cube): 4,288 problems. 160 and 384 are the first tilings that
  need more than one cube per cell (P = 2d k^(d-1) cubes at s = k, k = 2), so the first two containers have the tilings for k = 1 and 2 inside their ranges (4B.5). The other
  four stop at 400, below their second tiling (P = 896, 2,048, 4,608 and 10,240; the last is also above the id limit 9,999): a cap of 400 keeps the server's check under 20 ms and the kernel's quadratic pair stage small.
  The larger tilings are still certified by the tests. All 4,288 problems have lower bounds and no records.

### 4B.2 Lattice words in d - 1 axes

Verbatim from 4.1, with Z^3 replaced by Z^n. Every ridge is a right angle, so the cells a straight path crosses unfold into the cubes of the lattice 2h Z^n with
the anchor at the origin; the hinges are the lattice points with at least two coordinates that are odd multiples of h (their codimension in R^n is 2, as
a (d-3)-face's in the boundary). A **word** is the ordered list of the planes y_k = (2m + 1) h that a straight path crosses; it has an end cell and an exact
map x = A y + h b into R^d, A a signed permutation (stored per lattice axis k as a row r_k and a sign g_k, so that A e_k = g_k e_(r_k)) and b an integer vector.
Crossing in direction tau a plane of axis k, from cell (i, sg) with r_k = j and t = tau g_k, enters cell (j, t) by the rotation R of the plane (e_i, e_j) about the
ridge, R(sg e_i) = t e_j, R(t e_j) = -sg e_i; then A <- R A, b <- R (b - q) + q, q = sg e_i + t e_j.

**Completeness for the reach used** [CALC]. *Every straight path of length L < 2s is realised by a word with at most two crossings per lattice axis, all in one
direction.* Proof: three crossings of one axis are three consecutive planes, which are 4h = 2s apart in that coordinate (as in 4.1; the argument never used d).
Balls need paths shorter than 1 (their diameter) and fit only from s = 1/sqrt 2 (4B.4), so 2s > 1; from s >= 1 only one crossing per axis is needed (the third plane level
is at distance >= s from any centre). So the words with |m_k| <= 2 hold every path that can matter. Cubes need only one-crossing words (4B.5).

| d | n = d-1 | words per cell, <= 2 per axis | <= 1 per axis |
|---|---|---|---|
| 5 | 4 | 109,601 | 633 |
| 6 | 5 | 9,864,101 | 6,331 |
| 7 | 6 | 1,302,061,345 | 75,973 |
| 8 | 7 | 236,975,164,805 | 1,063,623 |
| 9 | 8 | 56,874,039,553,217 | 17,017,969 |
| 10 | 9 | 17,403,456,103,284,420 | 306,323,443 |

[CALC] The counts are sum_k (2n)!/(2n-k)! (two crossings: the exponential generating function (1+x)^(2n)) and sum_k C(n,k) 2^k k! (one). The first two rows are the research
note's 109,601, 9,864,101, 633 and 6,331; the tests enumerate d = 5 in full (109,601 words from three start cells), d = 6 in full (9,864,101, once) and the rest by formula. **The words are not listed**:
they are found for one pair of points at a time, by reach.

### 4B.3 The search by reach: the bounds, the loop lemma

For a start u (cell c1) and a target point x2 (cell c2) the search looks for the realised words of length < 1. It first chooses the **support** m (the crossings per axis, m_k in
{-2..2}, a path of length < 1 crosses axis k only up to the plane at distance < 1 from u), then the **orders** of the crossings, and tests each ordered word exactly as 4.1: the image
y = A^T (x2 - h b) lies in the closed lattice cube m, and the segment u -> y crosses the word's planes at parameters 0 <= t_1 <= ... <= 1 (a plane holding both ends puts no constraint).
Four exact facts prune it; each holds for every word of the support at hand, so the search finds every realised word of length < 1 (the Python tests compare it with a full enumeration, exact, d = 5 to 7).

1. **Support bound.** length^2 >= sum over axes with m_k != 0 of delta_k^2 + sum over axes with m_k = 0 of (x2[r_k] - u_k)^2, delta_k = (2|m_k| - 1) h - sign(m_k) u_k, the distance from u to the last plane
   crossed on axis k. The first sum: the segment crosses that plane, so |y_k - u_k| >= delta_k. The second: a column of A that is never crossed keeps its row r_k, and b keeps its value there, so y_k = x2[r_k].
2. **Assignment bound.** y = 2h m + yhat, where yhat is a signed permutation of the target's local coordinates (a column of A is a signed unit vector), so with c = 2h m - u,
   length^2 = sum (c_r + yhat_r)^2 >= sum (|c|_(r) - |u2|_(r))^2, both lists sorted ascending (the rearrangement inequality).
3. **Image coordinates are known early.** *A crossing of axis k changes only column k of A and the entries of b along the old normal axis i and the old row j.* (R fixes e_l for l != i, j; the columns of A are
   tangent to the cell, so only column k has a component along e_j, and none along e_i.) *No later crossing of another axis touches them:* a later crossing changes b along the then-normal axis and along a row
   of another column, and neither is the row that column k now holds (the normal and the rows are all distinct). Hence y_k = g_k (x2[r_k] - h b[r_k]) is fixed from the last crossing of axis k, and during the
   orders the search replaces delta_k by y_k - u_k in the sum, requires y_k in its lattice cube, and requires the crossing times known so far to lie in [0, 1] and to increase along the order.
4. **The end cell.** *After any word the normal axis is the original row of a crossed axis, or the start's own normal axis* (a crossing swaps the normal with the row of the crossed column; induction). So unless an end cell has the
   start's normal axis, the support must cross the axis whose original row is that end cell's normal axis.

The search is capped by a **work budget** (steps of the support and order searches, 4,000,000 per check, `NCUBE_MAX_WORK`; 8,000 per pair in the kernel's energy): a search over its budget throws, and the
check refuses the configuration with a plain reason. The cap makes the validator **sound but incomplete**: it never accepts an invalid configuration, and refuses a valid one only if
the balls are so crowded against the (d-3)-faces that the search needs more than its budget. [CALC] The worst pair found (two balls at the minimum clearance 1/2 near a common corner, every slack 0.3536 on all n axes) cost, measured and approximate, about 129, 374 and 1,105 steps for d = 5 to 7 and 3,294, 9,601 and 26,470 for d = 8 to 10 (an earlier figure for d = 8 to 10, 3,230, 9,025 and 23,590, undershot; none of these is an upper bound); the cap holds about 150 such pairs at d = 10 and takes about 0.25 s to exhaust.

**The loop lemma** [CALC]. *A nonempty word that ends in its start cell crosses at least two lattice axes twice.* Proof: write the frame (normal, A e_1, ..., A e_n) as the matrix F of a signed permutation;
crossing axis k in direction tau is F <- F G with G the 90 degree rotation of the plane (0, k), by tau. Conjugating by the diagonal matrix of the directions tau_k (which fixes e_0), every crossing becomes the
same rotation R_k: R_k e_0 = e_k, R_k e_k = -e_0, R_k e_l = e_l. The word returns to its start cell exactly when the product fixes e_0. Following e_0 through the product from the last crossing: it moves to e_k at a
crossing of k, stays until the next crossing of the same k, where it becomes -e_0 (the sign flips); and so on. Returning to +e_0 needs an even number of such excursions, two at least, each using
both crossings of one axis, and an axis crosses at most twice, so two different axes are crossed twice. [CALC: checked exhaustively for d = 5 and d = 6, all 109,601 and 9,864,101 words per cell.] So a geodesic loop is at least sqrt((3h - u_k)^2 + (3h - u_l)^2) >= sqrt 2 s
long (the last plane of an axis crossed twice is 3h - u_k >= 2h = s away), which is >= 1 once a ball fits, and the loop check can never fire for s >= 1/sqrt 2 (the validator and the certifier keep it anyway; the kernel's energy leaves it out).

### 4B.4 Balls

The rule is 4.2 unchanged (its proof is dimension-free): a packing of balls of diameter 1 is valid exactly when every centre is >= 1/2 from the hinges, every two centres are >= 1 apart along every realised
word, and every loop is >= 1.

- **Clearance**, closed form: min over pairs of axes a < b of sqrt((h - |u_a|)^2 + (h - |u_b|)^2), the root of the two smallest squared slacks (the nearest lattice hinge is an edge-of-cell hinge of the anchor
  cube; the Python tests check it against the lattice points). It is at most h sqrt 2, so a ball fits only from s = 1/sqrt 2 [CALC], and exactly there at the cell's centre.
- **Pairs.** Only pairs whose centres are < 1 apart in R^d are searched (the geodesic is never shorter than the chord), found in a grid on four coordinates, so the work is linear in P; at most 20 reasons.
- **No ridge snap.** Every gluing is an exact signed permutation, so no centre needs to be moved onto a facet (the chain validator of the 4-polytopes does). A centre on a ridge belongs to two cells: the search
  accepts either as the target's end cell.
- **Tolerance.** tol 1e-9 as the tesseract's (centres within tol of their cell, clearance >= 1/2 - tol, distances >= 1 - tol); the site's strict check for a problem without a record uses 1e-11.
- **Cost** [CALC, Apple Silicon, one core]: the engine checked layouts of one to twenty balls per cell, tight at their edges, in 0.8, 0.9, 1.4, 2.8, 6.3 and 9.0 ms (d = 5 to 10, P = 2d) and 5.4, 12.9, 10.5, 10.3,
  11.6 and 13.3 ms at P = 160, 384, 400, 400, 400, 400 (the test packings of the problem sets' largest P). The exact Python certifier takes 0.01 to 0.8 s for the P = 2d packings (d = 5 to 10).

### 4B.5 Axis-aligned cubes

4.3 with n axes. The cube centred at u is the box u +- 1/2 and crosses axis k when |u_k| + 1/2 > h. [CALC] *It is valid alone exactly when it crosses at most one axis* (two crossed axes k, l put the anchor's own hinge
u_k = +-h, u_l = +-h, the other coordinates u_m, inside the open box; for s < 1 every axis is crossed, so a cube needs s >= 1 and, with one crossed axis, the slab beyond the face has depth <= 1/2 and stays in the
neighbouring lattice cube). Its parts are the box clipped to the anchor and, beyond a crossed face, the slab mapped by the one-crossing word, a signed permutation, so again a box; parts of distinct cubes in a cell must
not overlap by more than tol in all n axes. Strict mode counts a face touched from inside as crossed. A cube never overlaps itself (its two parts lie in different cells). The word search is not needed.

**Tilings** [CALC]: P = 2d k^(d-1) cubes at s = k, k^(d-1) per cell on the grid of spacing 1, touching, so valid but not strictly valid; they are optimal for the cube problem by the volume bound (P <= 2d s^(d-1)),
which they meet, so no packing of cubes of any orientation does better [CALC]. By the policy for proofs they stay best_known unless a published proof exists (none is cited; no records ship).

### 4B.6 Lower bounds

Stored as the largest binary64 value not above the exact one, decided in integer arithmetic (`ncubeBound`), and checked exactly in Python (`ncube.py exact_bounds`): the larger of
- **volume**: P vol(piece) <= 2d s^(d-1) (the developments are injective on interiors, the interiors are disjoint); a ball has volume V_n 2^-n, V_n the unit n-ball's (pi replaced by its 60-digit truncation: a smaller, valid bound);
- **one piece**: 1/sqrt 2 for a ball (4B.4), 1 for a cube (4B.5).
For example P = 160 cubes on the 5-cube: exactly 2 (the tiling), and one ball: 0.7071067811865475 on every container.

### 4B.7 The exact certifier

`python/boundary/ncube.py`, reached through `certify.py` (`boundary/<d>cube/...` ids). Exact rationals only; shares no code with the engine. The same words, but a simpler search: the plane-distance part of bound 1 and bound 2, not the known image coordinates of fact 3 nor the end-cell fact 4 (and an unpruned mode, `prune=False`, that tries every word with
|m_k| <= 2): so it also checks the engine's pruning. Strict validity needs > in every test. Checks:
- the words (counts, every development an unfolding, the loop lemma), the pruned search against a full enumeration (exact, random pairs, d = 5 with two crossings and s from 0.75; d = 6 and 7 with one), the clearance against the lattice points;
- 42 test packings (`fixtures/ncube_packings.json`, built by `build_ncube_packings.py`): per d one ball at the cell's centre at s = (1/sqrt 2)(1 + 1e-9), balls at the centres of the 2d cells, two per cell (4d balls, on a main diagonal of the inner cube),
  one cube, 2d cubes, four cubes centred on the ridges of a ring of cells (each in two cells), and eight cubes round it, half a step off; each scaled by 1 + 1e-9 (4.4's argument) and certified strictly valid; their copies shrunk by 1e-6 fail;
- the exact tilings (d = 5, 6, 7, 8, 10, k = 1 and 2: up to 896 cubes at s = 2 exactly), 160 and 384 balls on the 2^n grids of the cells at s = 2, touching;
- the adversarial cases of both validators (`test_ncube.py`, `packages/boundary/test/ncube.test.ts`), d = 5 to 10: a ball over a (d-3)-face (clearance 0.509 valid, 0.495 not, on the hinge not), two balls across a ridge 1e-9 too close (the tolerance decides) and 1e-9 clear, two balls
  that overlap only through the word that crosses two ridges round a hinge (clearance 0.6 each, 85 and 185 degrees round it, path 0.919, the other way 1.195), a cube over a hinge, two cubes whose parts overlap across a ridge by 0 and 1e-9, malformed input;
- the engine's verdicts on 648 configurations (`fixtures/ncube_crosscheck.json`, `gen-ncube-crosscheck.ts`; the test packings scaled by 1.2, 1, 0.9999 and 0.97, and seeded random layouts, 325 valid): equal unless an exact margin lies within 1e-11 of its threshold.

### 4B.8 Server, kernel, caps

- **Rule** (`ncubeServerRule`): P at most the problem set's (160, 384, 400), the form (cell 0 to 2d - 1, u of n numbers at most 1e4 in absolute value, no `q`), s >= the bound less 1e-9 relative, then the validator at 1e-9 with 20 reasons.
  The space's P limit rises from 64 to 400 for this; each container's own range is its domain's.
- **Heavy checks.** `boundaryIsHeavy` uses 0.3, 0.35, 0.5, 0.6, 1 and 1.5 ms per ball (d = 5 to 10; three times the measured per-ball cost) and 0.02 ms per cube: a check of 50 ms or more waits in the server's heavy queue (one at a time) and is charged its estimate.
- **Caps.** At most 400 pieces, the work budget above (a worst-case check takes about 0.25 s), 20 reasons, cell and coordinate ranges, no `q`; a configuration that exhausts the budget is refused, never slow.
- **Kernel** (section 14): `boundaryKernel` returns `ncubeKernel` for these containers: n variables per piece, the exact energy of 4B.4 and 4B.5 with its analytic gradient in the displacements and in the size (the walk is `move.ts`, a signed permutation per crossing), weights, focus, a pair whose search exceeds
  8,000 steps in the energy falls back to its chord in R^d; seeds are the cell layouts (balls: one per cell centred at s = 1, two on a diagonal at s = 1 + 1/sqrt n, then the grids of spacing 1 at s = k; cubes: the k^n grids) and the **ring cut** (the cells of two axes form a
  flat strip of length 4s round the 2-axes' square, times a g^(d-2) grid across it, pieces at equal spacing from a random phase, some straddling ridges; a ball ring is kept only if the validator accepts it); finishing is the search's bisection plus the
  validator at tolerance 0. The generic methods run on it; `searchBoundary` on a d-cube problem seeds and finishes (no optimiser: no records are searched).
- **Tighter** must not see these containers: the game's slices and id patterns do not name `boundary` (`apps/web/test/tighter-boundary-excluded.test.ts` pins that for every d-cube id).

### 4B.9 Limitations

- **No records.** Test packings only. The cell layouts and the tilings are exact or nearly so; the others have no claim to being good (no search ran, as the tesseract's first two rounds would be); the optima for balls are open; a tesseract-like P = 2d layout below s = 1 would need the pair search at the worst corner (cost above).
- **Aligned cubes only**, as on the tesseract: a turned cube might do better and would need parts and exact rotations in n dimensions (`parts.ts`, `polyhedron.ts` are 3-dimensional).
- **The corner cap.** The validator refuses a configuration whose search needs more than 4,000,000 steps; the worst pair of the kind above costs about 26,000 at d = 10 (measured, not a bound), but a few hundred balls crowded at a common corner of d = 10 could reach it.
- **The vertex rule** (section 3) is argued, not proved, in every d.
- Not covered: the simplex, orthoplex and demicube boundaries for d >= 5, whose dihedral angles are not right angles; the chain code would need to go to R^(d-1) (phase C3 of the plan).


## 4C. Platonic solids on the regular 4-polytopes

Added 2026-11-02. Pieces other than balls: the tetrahedron, cube, octahedron, dodecahedron and icosahedron of **edge 1**, freely rotated, on the
boundary of the tesseract, 5-cell, 16-cell, 24-cell, 120-cell and 600-cell (4A.1). Same score (the container's edge s), same gluing,
same rule (sections 2 and 3: no edge of the polytope inside a developed piece, the pieces' interiors pairwise disjoint, a piece not overlapping
itself), same tolerances as the balls and the aligned cubes. This section adds what the solids need: their coordinates, the validator's parts, the
exact certifier (a field of degree 4, an exact rotation), the lower bounds, the search with rotations and the test packings.

### 4C.1 The pieces, the ids, the problems

- **Standard solids** (`packages/boundary/src/pieces.ts`; `python/boundary/platonic.py` repeats them exactly), unit edge, centred at the origin, phi = (1 + sqrt 5) / 2:
  the tetrahedron (1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1) / (2 sqrt 2); the cube (+-1/2)^3; the octahedron +-e_k / sqrt 2; the dodecahedron (+-phi/2)^3 and the
  cyclic shifts of (0, +-1/2, +-phi^2/2); the icosahedron the cyclic shifts of (0, +-1/2, +-phi/2). Every Platonic solid is mirror-symmetric, so a
  rotation is all a piece needs (an improper orientation, which a gluing that reverses handedness produces, is the same piece turned by the solid's mirror,
  `PieceShape.mirror`).
- **Pose.** `{cell, u, q}` as for the cubes (section 2): q is the unit quaternion (w, x, y, z) that rotates the standard solid about its centre in the anchor's
  local frame; absent means no rotation. The circumradius (0.612, 0.866, 0.707, 1.401, 0.951 for the five) is how far a piece reaches from its centre.
- **Ids.** `boundary/<container>/<piece>/<P>` with piece `tetrahedron`, `cube`, `octahedron`, `dodecahedron` or `icosahedron`. **The id `cube` on the 5-cell, 16-cell, 24-cell,
  120-cell and 600-cell is a freely rotated cube; on the tesseract it is the axis-aligned cube of section 4.3, with its records, rules and ids unchanged.**
  The tesseract takes the four other solids freely rotated. The d-cubes (4B) take balls and aligned cubes only.
- **P** = 1 to 24 on the five general containers, 1 to 16 on the tesseract, as for their balls.
- **Problems.** 664 more than before: the tesseract's four solids (4 x 16 = 64) and five solids on each of the five others (5 x 5 x 24 = 600); 816 on the 4-polytopes in all (152 + 664),
  none with a research record (status `no_record`). A first entry on such a problem must pass the strict tolerance 1e-11 of section 12.

### 4C.2 The float validator

`validateBoundary` (`packages/boundary/src/validate.ts`) is the cube validator made general. A solid X with pose (u, R) is developed from its centre exactly as a cube: the chains of cells
within its circumradius are enumerated (chains.ts), and for each chain C ending in cell c the part of X in c is X cut by c's facets and by C's window, mapped
into c's frame (`partPolyhedron` clips the solid's faces; the piece's face normals are carried through the placement, so a reversed handedness stays outward). The checks:

1. **Edges of the polytope.** Every hinge of every chain's cell, clipped to the closed window, against the interior of the solid: the depth of a point x in the solid is the least over the
   faces of (inradius - n_k . (R^T (x - u))), concave and piecewise linear along a segment, so its maximum is at an end or where a rising face function crosses a falling one
   (`segmentDepthInShape`); the cube keeps its closed form. A depth above the tolerance rejects.
2. **Parts.** In every cell all parts, those of one piece included, have disjoint interiors: a uniform grid of side max(2, diameter) per cell bins the parts' bounding boxes, and two
   parts whose boxes overlap by more than the tolerance get the separating-axis test on their face normals and edge-edge axes (`satDepth`); its depth is the penetration.
3. **The volume bound** first: P V_piece <= V_boundary s^3 (1 + 1e-6), then the in-cell edge tests, so that a container far too small for its pieces fails before any chain is enumerated.

*Work caps (fail closed).* At most 1024 chains per piece on the server (4096 for the library), a configuration that reaches it is rejected with "reaches more than N chains of cells; the
polytope is far too small for its pieces"; at most 400,000 tests of pairs of parts (`maxPairs`), a configuration that needs more is rejected with the reason "more than N tests of pairs of parts; the
configuration is too crowded for the server to check". A cap never turns a failure into a pass.

*Measured times* [CALC, 2026-11-02, Apple Silicon, one core, tolerance 1e-9, 1024 chains]: one solid alone at 1.02 to 4 times its smallest edge in 12 random orientations, at most 2.2 ms
(single outliers of 6 ms on the tesseract and 3.8 ms on the 120-cell); a loose random layout of P = 24 (16 on the tesseract) 1.4 to 25 ms; the cells themselves, every solid touching its
neighbours, 3.8 ms a piece on the 5-cell, 2.8 on the 16-cell, 2.3 on the 24-cell, 4.4 on the 120-cell (120 dodecahedra, 0.53 s) and 1.2 on the 600-cell (600 tetrahedra, 0.73 s).
**Heavy rule.** The server prices a solid at `BOUNDARY_SOLID_MS` per piece (the worst of these, rounded up: 4 ms on the tesseract and the 5-cell, 3 on the 16-cell and 24-cell, 5 on the 120-cell,
2 on the 600-cell; the aligned cubes keep 0.3 and the balls their own), and a check estimated at 50 ms or more waits in the heavy queue: solids from P = 13 on the tesseract and the 5-cell,
17 on the 16-cell and 24-cell, 10 on the 120-cell, never on the 600-cell (24 x 2 ms).

### 4C.3 The exact certifier

`python/boundary/platonic.py` (called by `certify.py` for these problems through the same entry point; `check_pieces`) reuses `general.py`'s container, charts, gluings and chain enumerator unchanged
and adds `python/boundary/qmulti.py`, exact arithmetic in Q(sqrt 2, sqrt 5) (elements (n0 + n1 sqrt 2 + n2 sqrt 5 + n3 sqrt 10) / d in integers; signs decided exactly, with a rigorous float
fast path; `Q5` accepts it in mixed arithmetic). It shares no code with the TypeScript engine.

*The field.* A piece's vertex is x + sum_k (kappa / s) (R v)_k f_k (R an exact rational rotation, v the standard vertex, f_k three mutually orthogonal vectors of equal squared length m
of the cell's hyperplane, built from the cell's own vertices: the tetrahedron's circumscribed-cube axes, the octahedron's diameters, the cube's edges, the dodecahedron's inscribed-cube edges), where kappa = e0 / sqrt m
is 1/sqrt 2 (tetrahedral and octahedral cells), 1 (cubes) or 1/phi (dodecahedra) and e0 the base polytope's edge (`cell_frame` checks kappa^2 m = e0^2 in the field). The coefficients' field,
and so the pieces' field beyond the container's own (Q for the tesseract, 5-cell, 16-cell, 24-cell; Q(sqrt 5) for the 120-cell and 600-cell):

| | tetrahedron | cube | octahedron | dodecahedron | icosahedron |
|---|---|---|---|---|---|
| tesseract | Q(sqrt 2) | (aligned: 4.3) | Q(sqrt 2) | Q(sqrt 5) | Q(sqrt 5) |
| 5-, 16-, 24-cell | Q | Q(sqrt 2) | Q | Q(sqrt 2, sqrt 5) | Q(sqrt 2, sqrt 5) |
| 120-cell | Q(sqrt 2, sqrt 5) | Q(sqrt 5) | Q(sqrt 2, sqrt 5) | Q(sqrt 5) | Q(sqrt 5) |
| 600-cell | Q(sqrt 5) | Q(sqrt 2, sqrt 5) | Q(sqrt 5) | Q(sqrt 2, sqrt 5) | Q(sqrt 2, sqrt 5) |

The cells' own frames matter. The plan of 4A.8 (written before this) took the tetrahedron and octahedron to need the degree-4 field on both the 120-cell and the 600-cell. It holds for the 120-cell (the frame of a
dodecahedral cell has squared length phi^2 e0^2, so 1/sqrt 2 stays in the coefficients). The 600-cell's tetrahedral cells have the circumscribed-cube frame of squared length 2 e0^2, and kappa = 1/sqrt 2 cancels the 1/(2 sqrt 2) of the
tetrahedron's and the 1/sqrt 2 of the octahedron's coordinates: Q(sqrt 5) there, and Q on the 5-cell, 16-cell and 24-cell. [CALC: `Frames` in `test_platonic.py` checks every frame exactly; the table is read off the vertices.]

*Reading a configuration* (the engine's floats). The centre is read as in 4A.3 (dyadic weights of four vertices, an exact point of the cell's hyperplane, moved onto the facet planes it exceeds by less than
1e-12 of an edge). The rotation: the engine's frame F_e and the exact frame F_x (scaled to unit length, in the engine's float coordinates by the affine map that takes the cell's exact vertices to the engine's
unit-edge vertices) span the same 3-space, so the engine's rotation in the exact frame is O R_e, O = F_x^T F_e, R_e the matrix of the normalised float quaternion; if det O < 0 the third exact frame vector changes sign.
That float matrix becomes a quaternion, rounded to the grid 2^-60, and the exact rational rotation of that rational quaternion (entries rational, exactly orthogonal, determinant 1) is the pose's. A centre that lies on
a ridge is developed from every closed cell that contains it (the runs' parts are compared together), as for the balls.

*Discrepancy and slack.* [CALC] Over 540 random poses (every solid on every container, s = 0.6, 2 and 5, centres anywhere in their cells, random rotations) and the 12 test packings, the largest distance between a
vertex of the exact pose and the engine's own float computation of it is **2.74e-15 s** (120-cell cube), 1.2 to 2.3e-15 s on the others; for balls 4A.3 found 2.9e-15 s. Strict validity is decided on the
pieces **inflated about their centres by the factor 1 + EPS_STRICT, EPS_STRICT = 2^-40 = 9.1e-13**: if the inflated solids are valid (interiors disjoint, no edge through an interior, embedded), the closed solids are
strictly valid, each at distance at least EPS_STRICT times its inradius from the others and from every edge, which is 1.9e-13 for the tetrahedron (the smallest) against a reading error of 2.74e-15 s: a factor 70 at s = 1 and 7 at
s = 10. The certifier refuses (with that reason, never accepts) a strict check with s READING_BOUND > EPS_STRICT r_in, READING_BOUND = 1e-14 (3.6 times the measured error): s above 18 for the tetrahedron, 37 for the
octahedron, 45 for the cube, 68 for the icosahedron, 101 for the dodecahedron. The reach of the chain enumeration is the inflated circumradius.

*Decisions* (exact): the lower bounds first (4C.4; strict: s equal to the bound is invalid); the centre in its closed cell; then for each piece and run the chains within the circumradius (at most 40,000 per piece and run,
refused with a reason beyond), and in each node (1) every hinge of the node's cell clipped to the closed window against the open interior of the piece (an interval test per face), (2) the part: the piece's vertices
in the cell's own chart, clipped one plane at a time by the cell's facets and the window's planes (exactly, each vertex keeping the bitmask of the planes it lies on, so the part's edges are known without ordering
anything; parts without interior are dropped), and (3) in every cell every two parts of different (piece, run, chain), parts of one piece included, interior-disjoint: the bounding boxes first (rigorous floats, margin 1e-7),
then the separating-axis test, a pair being accepted **only with an exact separating plane** among the face normals of either and the cross products of an edge of each (candidates ordered by their float gap, the
first 12 face and 40 edge candidates verified exactly). A pair without a verified plane is reported as overlapping, so an unexpected case is refused, never accepted.

*Cost* [CALC, Apple Silicon, a loaded machine]: certifying 2 or 3 solids takes 0.1 to 2.6 s (2.6 s for 2 dodecahedra on the 120-cell); the 5 tetrahedra of the 5-cell 0.8 s, the 16 of the 16-cell 2.6 s, the 24 octahedra of the 24-cell 3.3 s, 24 tetrahedra of the 600-cell 4.8 s, 24 dodecahedra of the 120-cell 25 s.

### 4C.4 Lower bounds

Both bounds of section 5 hold for solids [CALC]; stored as binary64 rounded down by a relative 1e-15 (`pieceBounds`, `packages/boundary/src/bounds.ts`), checked exactly by `platonic.exact_bounds`
(enclosures of sqrt 2 and sqrt 5 to 45 digits; `Bounds` in `test_platonic.py` checks all 664 stored values: each not above the exact bound and within 1e-12 of it).

- **Volume.** The developed solids have disjoint interiors and the developments are injective, so P V_piece <= V s^3, V the volume of the boundary at unit edge (4A.4). V_piece = sqrt 2 / 12 (tetrahedron), 1 (cube),
  sqrt 2 / 3 (octahedron), (15 + 7 sqrt 5) / 4 (dodecahedron), 5 (3 + sqrt 5) / 12 (icosahedron).
- **One piece.** The solid contains the open ball of its inradius r about its centre (r = sqrt 6 / 12, 1/2, sqrt 6 / 6, sqrt (250 + 110 sqrt 5) / 20, phi^2 / (2 sqrt 3)), and the development of that ball is part of the solid's: it is injective and
  contains no point of an edge. By 4.2 the centre has clearance at least r from the edges; but every point of a cell lies within s rho of the cell's edges, rho the cell's midradius (4A.4). So s >= r / rho. For the cube
  on a container other than the tesseract this is the ball bound of 4A.4 (r = 1/2).

One-piece bound / volume bound at P = 24 (16 on the tesseract), 6 decimals; the larger is the stored bound:

| container | tetrahedron | cube | octahedron | dodecahedron | icosahedron |
|---|---|---|---|---|---|
| tesseract | 0.288675 / 0.617715 | 0.707107 / 1.000000 (aligned, 4.3) | 0.577350 / 0.980561 | 1.574750 / 2.483963 | 1.068808 / 1.634085 |
| 5-cell | 0.577350 / 1.686865 | 1.414214 / 3.440613 | 1.154701 / 2.677732 | 3.149500 / 6.783248 | 2.137616 / 4.462386 |
| 16-cell | 0.577350 / 1.144714 | 1.414214 / 2.334815 | 1.154701 / 1.817121 | 3.149500 / 4.603142 | 2.137616 / 3.028195 |
| 24-cell | 0.408248 / 0.629961 | 1.000000 / 1.284898 | 0.816497 / 1.000000 | 2.227033 / 2.533207 | 1.511523 / 1.666480 |
| 120-cell | 0.155937 / 0.145430 | 0.381966 / 0.296625 | 0.311874 / 0.230855 | 0.850651 / 0.584804 | 0.577350 / 0.384715 |
| 600-cell | 0.577350 / 0.341995 | 1.414214 / 0.697550 | 1.154701 / 0.542884 | 3.149500 / 1.375236 | 2.137616 / 0.904705 |

*Where a bound is met* [CALC]: the 5 tetrahedra of the 5-cell, the 16 of the 16-cell, the 24 octahedra of the 24-cell and (as a subset of the cells) any number of dodecahedra of the 120-cell and tetrahedra of the
600-cell fill their cells at s = 1; for P = 5, 16 and 24 the volume bound is exactly 1 on the first three, so those packings are optimal by our own derivation [CALC]. PolyPack calls a value proven only with a
published proof (section 5), so they stay unrecorded here: these problems have no records at all. Beyond them nothing is known about optimality [SPEC].

### 4C.5 Checks around the certifier

- **Engine against certifier** (`fixtures/platonic_crosscheck.json`, `gen-platonic-crosscheck.ts`): for every container and solid (29 pairs), two solids at random cells and turns, a third of them with a centre on a ridge, a
  second family at the centres of two adjacent cells; the engine's validator finds the critical scale k* of the whole configuration (bisection on validity, monotone in the scale), and the configuration is emitted at k* (1 +- d), d = 1e-5
  everywhere and 1e-9 for the tolerance-0 entries of every fourth pair. 58 entries, 136 configurations: 11 with a centre on a ridge, 29 at adjacent centres; the refused side's reason is an edge through a solid in 44 and
  an overlap in 24. **The exact verdict (plain validity) is the engine's on all 136** (tolerance 0, or 1e-9 for the entries with a centre on a ridge, whose float coordinates can be a rounding unit outside their cell), and the
  strict verdict on the 1e-5 side of every fourth entry agrees too. So the float validator never accepts what the exact certifier refuses, and the reverse, outside a margin of 1e-9 of s.
- **Frames and field.** Every frame is orthogonal, of equal length m, with kappa^2 m = e0^2; its image in the engine's coordinates is orthonormal and orthogonal to the engine's normal (every cell of the 5-, 16- and 24-cell and the tesseract,
  six of the 120-cell and 600-cell); the field's signs against enclosures; the rational rotations exactly orthogonal; the tesseract's cubes by this path against the lattice words of 4.3 (s = 1 +- 1e-6 and exactly 1, plain and strict).
- **Adversarial tests** (`test_platonic.py`, `packages/boundary/test/solids.test.ts`): a tetrahedron filling its cell pushed through a face; the octahedra of the 24-cell at s = 1 + 1e-6, 1 + 1e-9 (valid, strictly valid) and 1 - 1e-6 (an edge through
  them); touching variants (the cells at s = 1 are valid, not strictly); a centre on a ridge, and on an edge (refused); malformed input; the work caps; an s beyond the strict reading bound; the volume bound.
  A solid cannot overlap itself around a sharp edge: a convex piece spans at most 180 degrees round an edge, and the smallest total angle is 211.6 (5-cell); 4000 random poses on each of four containers and five solids found no
  such overlap, only edges through a solid.

### 4C.6 The search, the kernel and the test packings

*Model* (`packages/boundary/src/search/solids.ts`, `SolidModel`). A state is each piece's anchor, centre u and orientation R (3 x 3, orthogonal, improper allowed). The penalty energy is the sum of squared
penetration depths: (1) for each pair of pieces and each chain of cells from i's centre within twice the circumradius whose window holds the image of j's centre, the separating-axis depth of the rigid copy of j
(centre T u_j, orientation M R_j) against i (a copy of a piece along a non-trivial chain is a loop and counts alike); (2) for each edge of the polytope, its deepest hit inside a piece over the chains. It is the validator's rule with
rigid copies instead of clipped parts, equal where the development is flat; the final scale and nudge (`finalizeState`) are decided by the validator, so a difference costs quality, never validity. Gradients are analytic in the
centres, in the rotation increments theta (R' = exp(theta) R, theta in the anchor's frame) and in s at fixed centres, by the contact feature that sets the depth (a face against a vertex, an edge against an edge, the active
face(s) of a hinge hit). [CALC] They agree with central differences on all six containers and five solids (to 1e-9 in the checks run; `solids-search.test.ts` asks 1e-6). Cost of an energy evaluation with its gradient per piece, the
worst of the five solids on squeezed random starts: 0.2 to 0.9 ms on the tesseract, 5-cell, 16-cell, 24-cell and 600-cell, 3.4 ms on the 120-cell (`BOUNDARY_SOLID_ENERGY_MS` rounds these up: 0.7, 0.4, 0.35, 1, 0.3 and 4).

*Driver.* `searchBoundary` adds three rotation variables per piece after the 3P displacements: a displacement is moved through the charts, the orientation after it is M exp(omega) R, and the gradient in omega is
Jl(omega)^T det(M) M^T times the model's gradient in theta (Jl the left Jacobian of the exponential, det M = -1 when a gluing reverses handedness). Starts: the cells themselves when the solid is the container's cell and P <= C
(`constructedStart`, valid from s = 1), random layouts at a start edge at which a cell's inscribed ball holds the solid's circumscribed ball; hops move one to three pieces.

*Kernel* (`packages/boundary/src/kernel-solids.ts`, section 14): `dim` = 6 (the displacement and a rotation vector about the anchor orientation), a typical step the inradius and inradius / circumradius radians, `swap`
exchanges two pieces' orientations (an orientation of the wrong handedness for the walk is replaced by its mirror image, the same piece), `seeds` as above, `finishCost` 60 and `validateCost` 15 energy evaluations. The
generic methods' contract (`checkKernel`) passes on every container and solid tested.

*Test packings* (`python/boundary/fixtures/platonic_packings.json`, 12; raw configurations from `gen-platonic-packings.ts`, scaled by lam_k = 1 + k 1e-12 about the cells' centres to the first k at which the exact
certifier accepts them strictly **and with twice the inflation**, then stored with the largest power of two inflation they pass): the cells themselves (5 tetrahedra of the 5-cell, 16 of the 16-cell, 24 octahedra of the 24-cell,
24 of the 120-cell's dodecahedra, 24 of the 600-cell's tetrahedra) and one short seeded run per container and solid. They are not good packings (no long search ran) and no optimality is claimed.

| problem | source | s (score_cert) | lower bound | k | inflation certified (margin) | reading error | chains / parts |
|---|---|---|---|---|---|---|---|
| 5cell/tetrahedron/5 | the cells themselves, s = 1 | 1.000000000002 | 1.000000 (volume) | 2 | 9.1e-13 (1.8e-12) | 3.0e-16 | 385 / 5 |
| 16cell/tetrahedron/16 | the cells themselves, s = 1 | 1.000000000002 | 1.000000 (volume) | 2 | 9.1e-13 (1.8e-12) | 5.4e-16 | 1232 / 16 |
| 24cell/octahedron/24 | the cells themselves, s = 1 | 1.000000000002 | 1.000000 (volume) | 2 | 9.1e-13 (1.8e-12) | 4.2e-16 | 792 / 24 |
| 120cell/dodecahedron/24 | 24 of the 120 cells, s = 1 | 1.000000000002 | 0.850651 (one piece) | 2 | 9.1e-13 (1.8e-12) | 1.7e-15 | 1752 / 24 |
| 600cell/tetrahedron/24 | 24 of the 600 cells, s = 1 | 1.000000000002 | 0.577350 (one piece) | 2 | 9.1e-13 (1.8e-12) | 4.5e-16 | 1848 / 24 |
| tesseract/icosahedron/2 | search, seed 1, 40000 | 2.137615801704 | 1.068808 (one piece) | 2 | 9.1e-13 (1.8e-12) | 3.3e-16 | 11 / 7 |
| tesseract/tetrahedron/3 | search, seed 2, 90000 | 0.577288575347 | 0.353553 (volume) | 2 | 9.1e-13 (1.8e-12) | 4.4e-16 | 204 / 18 |
| 5cell/octahedron/2 | search, seed 3, 60000 | 1.561196353847 | 1.169607 (volume) | 2 | 9.1e-13 (1.8e-12) | 2.1e-16 | 41 / 12 |
| 16cell/cube/2 | search, seed 4, 60000 | 2.103033866188 | 1.414214 (one piece) | 0 | 9.1e-13 (5.8e-11) | 2.3e-16 | 46 / 11 |
| 24cell/tetrahedron/3 | search, seed 5, 60000 | 0.800000000033 | 0.408248 (one piece) | 2 | 9.1e-13 (1.8e-12) | 4.2e-16 | 387 / 15 |
| 120cell/dodecahedron/2 | search, seed 6, 40000 | 1.000000000002 | 0.850651 (one piece) | 2 | 9.1e-13 (1.8e-12) | 1.7e-15 | 146 / 2 |
| 600cell/octahedron/2 | search, seed 7, 60000 | 1.250000000224 | 1.154701 (one piece) | 2 | 9.1e-13 (1.8e-12) | 1.3e-15 | 52 / 13 |

### 4C.7 Limitations

- No research records: the 664 new problems have lower bounds and the test packings above only. Nothing is known about their optimal values [SPEC]; the test packings are not good ones (no long run) and nothing about their
  distance from the optimum is claimed.
- The model (4C.6) is not the validator: where a piece wraps round an edge, rigid copies and clipped parts differ, so a search can be slower or stop higher than the exact rule allows; it cannot return an invalid packing.
- The search is slow and its first results are rough: the first valid packing (a random start relaxed and drawn in to touching) comes after 700 evaluations for one solid and 5,000 for six on the tesseract, but the compression
  that follows needs tens of thousands (40,000 to 90,000 for the test packings), at 0.1 to 4 ms per piece and evaluation. The board sizes its first round accordingly (section 13).
- The exact certifier is slow on crowded packings (25 s for 24 dodecahedra of the 120-cell; 40,000 chains per piece and run is its cap).
- The vertex rule (section 3) is argued, not proved, on every container, as before.
- A configuration with a centre exactly on an edge is refused as invalid (the edge is inside the piece).
- Whether the tesseract should also have a rotated-cube problem (its `cube` is the aligned cube) is a question for the project owner; the id `boundary/tesseract/cube/P` is the aligned problem and stays so.

## 5. Lower bounds

Only proven bounds are stored [CALC]. (The Platonic solids on the 4-polytopes have the two bounds of 4C.4.) Each is stored as a binary64 value rounded down, and the
record's `bound` is the larger one.

- **Volume.** Developments are local isometries, injective on interiors, with disjoint images, so
  the pieces' volumes add up to at most the boundary's, 8 s^3. Balls: s >= (P pi / 48)^(1/3).
  Cubes: s >= (P / 8)^(1/3).
- **One piece.** Every point of a cell lies within s / sqrt 2 of the cell's edges (the clearance
  formula is at most h sqrt 2 = s / sqrt 2, reached at the centre). A ball's centre needs
  clearance >= 1/2, so s >= 1/sqrt 2. A cube contains the ball of radius 1/2 about its centre, which
  must avoid the edges, so s >= 1 / sqrt 2.

Balls P = 1 and 2 meet 1/sqrt 2, and cubes P = 8 meet the volume bound 1 (the 8 cubes fill the
boundary). These values are therefore optimal by our own derivation [CALC]. PolyPack calls a value
proven only with a published proof, so they stay `best_known`, with a note. Every other value,
the search values and the exact constructions included, lies above both bounds; the tests check
this exactly.

"No unit cube fits on a tesseract of edge below 1" is numerical only (research note sec. 4:
maximising the side of a cube that avoids the edges, over rotation and centre, by 150 Nelder-Mead
runs and 200,000 random poses, gave 1) and is not used as a bound.

## 6. nofold values

[CALC, research note sec. 2.3] nofold(P) = f(ceil(P / 8)), with f(n) the smallest cube holding n
pieces. Upper bound: an optimal n-packing in every cell. Lower bound: some cell holds at least
ceil(P / 8) pieces inside it, which is an ordinary 3D packing.

| pieces | P = 1..8 (one per cell) | P = 9..16 (two per cell) |
|---|---|---|
| balls | 1 [CALC: a ball of diameter 1 in a cube needs edge 1] | 1 + 1/sqrt 3 = 1.577350 [CALC: the two centres lie in the inner cube of edge s - 1, of diameter sqrt 3 (s - 1) >= 1; equality on a diagonal] |
| cubes | 1 [CALC] | 2 [EST: two unit cubes in a cube, the cube-in-cube table of polyhedrainpolyhedra.wiki, as cited by the research note; upper bound CALC] |

## 7. Results

fold = score_cert at 6 decimals. Gain = 1 - score_float / nofold. Lower bound as in section 5.
Source: exact = a construction of the research note; research = the research prototype's round-2
search (2026-10-24; balls seeds 21 and 22, cubes seed 31 or 32); search = the PolyPack boundary
search (2026-10-25, one deterministic call with seed 1, below). "made exact" = section 4.6.

**Balls:**

| P | fold | nofold | gain | lower bound | source |
|---|---|---|---|---|---|
| 1 | 0.707107 | 1 | 29.3% | 0.707107 (one ball) | exact: centre of cell 0 |
| 2 | 0.707107 | 1 | 29.3% | 0.707107 (one ball) | exact: centres of cells 0 and 1 |
| 3 | 0.840071 | 1 | 16.0% | 0.707107 (one ball) | research (seeds agree) |
| 4 | 0.853553 | 1 | 14.6% | 0.707107 (one ball) | research, made exact: (2 + sqrt 2)/4 |
| 5 | 0.893613 | 1 | 10.6% | 0.707107 (one ball) | search |
| 6 | 0.932917 | 1 | 6.7% | 0.732296 (volume) | research (seeds agree) |
| 7 | 0.962702 | 1 | 3.7% | 0.770907 (volume) | research (seeds differ by 2.6e-4) |
| 8 | 0.965926 | 1 | 3.4% | 0.805996 (volume) | search, made exact: cos 15 deg |
| 9 | 1.110298 | 1.577350 | 29.6% | 0.838270 (volume) | research (seeds agree) |
| 10 | 1.134634 | 1.577350 | 28.1% | 0.868233 (volume) | search |
| 11 | 1.157959 | 1.577350 | 26.6% | 0.896260 (volume) | search |
| 12 | 1.164086 | 1.577350 | 26.2% | 0.922635 (volume) | research (seeds agree) |
| 13 | 1.198748 | 1.577350 | 24.0% | 0.947583 (volume) | search |
| 14 | 1.217715 | 1.577350 | 22.8% | 0.971282 (volume) | research (seeds agree) |
| 15 | 1.251620 | 1.577350 | 20.7% | 0.993879 (volume) | search |
| 16 | 1.263287 | 1.577350 | 19.9% | 1.015491 (volume) | search |

Mean gain 19.5%. "Seeds agree" means the research search's two independent seeds differ by at
most 1e-6 relative: evidence of the best basin found, not a proof. The gain falls from 29% to 3%
as P grows to 8, then jumps back at P = 9, where nofold must put two balls in one cell.

**Cubes** (axis-aligned):

| P | fold | nofold | gain | lower bound | source |
|---|---|---|---|---|---|
| 1..2 | 1.000000 | 1 | 0% | 0.707107 (one cube) | exact: one cube per cell |
| 3..7 | 1.000000 | 1 | 0% | 0.721125 .. 0.956466 (volume) | exact: one cube per cell |
| 8 | 1.000000 | 1 | 0% | 1 (volume) | exact: the cells themselves; meets the bound |
| 9 | 1.250000 | 2 | 37.5% | 1.040042 (volume) | research, made exact: 5/4 |
| 10 | 1.250000 | 2 | 37.5% | 1.077217 (volume) | P = 11 minus a cube: 5/4 |
| 11 | 1.250000 | 2 | 37.5% | 1.111990 (volume) | search, made exact: 5/4 |
| 12 | 1.250000 | 2 | 37.5% | 1.144714 (volume) | search, made exact: 5/4 |
| 13 | 1.333333 | 2 | 33.3% | 1.175667 (volume) | search, made exact: 4/3 |
| 14 | 1.333333 | 2 | 33.3% | 1.205071 (volume) | search, made exact: 4/3 |
| 15 | 1.500000 | 2 | 25.0% | 1.233106 (volume) | P = 16 minus a cube: 3/2 |
| 16 | 1.500000 | 2 | 25.0% | 1.259921 (volume) | research, made exact: 3/2 |

Mean gain 16.7%. Every cube value from P = 9 on is an exact construction at 5/4, 4/3 or 3/2
[CALC]; whether these are optimal is open [SPEC].

**What the search changed** (2026-10-25; `searchBoundary({polytope: "tesseract", piece, P},
{seed: 1, budget: {evals: E}})`, reproduced by `packages/boundary/scripts/search-sources.ts`).
It beat the research values at 11 problems. Improvement = 1 - new score_cert / research score_cert:

| problem | evaluations E | research | now | improvement |
|---|---|---|---|---|
| ball 5 | 2,617,672 | 0.893624 | 0.893613 | 0.0012% |
| ball 8 | 65,479 | 0.968315 | 0.965926 (cos 15 deg) | 0.25% |
| ball 10 | 654,968 | 1.135244 | 1.134634 | 0.054% |
| ball 11 | 1,358,135 | 1.157975 | 1.157959 | 0.0014% |
| ball 13 | 7,650,384 | 1.198812 | 1.198748 | 0.0054% |
| ball 15 | 5,026,907 | 1.252406 | 1.251620 | 0.063% |
| ball 16 | 779,083 | 1.264053 | 1.263287 | 0.061% |
| cube 11 | 156,372 | 1.333333 | 1.250000 (5/4) | 6.25% |
| cube 12 | 490,433 | 1.500000 | 1.250000 (5/4) | 16.7% |
| cube 13 | 298,374 | 1.500000 | 1.333333 (4/3) | 11.1% |
| cube 14 | 2,007,977 | 1.500000 | 1.333333 (4/3) | 11.1% |

The exact constructions also moved balls P = 4 (by 5e-11 relative, to (2 + sqrt 2)/4) and cubes
P = 9, 10, 15, 16 (by 1e-7 to 3e-7, onto 5/4 and 3/2). Each record keeps the research value in
its field `previous`.

## 8. Exact and float

- **Exact:** the certification of every record (strict validity at score_cert); the exact
  constructions of section 4.6 (rational, or in Q(sqrt 2, sqrt 3)); the lower bounds; the closed
  forms of the nofold values; the word maps.
- **Float:** the search values (score_float) of the ball records that are not exact
  constructions; the engine's validator; the cross-check numbers.
- The exact constructions sit on their exact values: balls P = 1, 2 at 1/sqrt 2, P = 4 at
  (2 + sqrt 2)/4, P = 8 at cos 15 deg; cubes P = 1..8 at 1, P = 9..12 at 5/4, 13, 14 at 4/3, 15, 16
  at 3/2. Their records store the binary64 rounding, scaled to strict validity, so score_cert is
  about 1e-12 above the exact value.

## 9. Limitations

- The vertex rule is argued, not proved [SPEC] (section 3). It does not affect the records.
- Ball values are float search results, except the exact constructions. The research value of
  P = 7 is not converged (two seeds differ by 5.1e-4), and the search values (P = 5, 10, 11, 13,
  15, 16) come from one seeded run each. All are upper bounds. Nothing is claimed optimal beyond
  the three cases that meet their bounds, and those stay `best_known`.
- Cubes: only axis-aligned cubes were searched. Rotated cubes were tried only by random
  single-cube insertion in the research (no fit found). The certifier refuses rotated cubes.
- The tesseract is the only container with research records, and the only one with cubes. The five others take balls (section 4A); their
  problems have lower bounds and test packings, and nothing better is known. Their Platonic pieces are not built (4A.8).
- The cube nofold value 2 for P = 9..16 is taken from the wiki's table, not re-derived here.
- The certifier's soundness rests on the arguments of section 4 and on Python's exact integers
  and fractions (and, for Q(sqrt 2, sqrt 3), on the sign rule of 4.6). Its records are also
  accepted by the engine's independent float validator.

## 10. Files and commands

- `python/boundary/certify.py`: the certifier (library and records mode);
  `python/boundary/qfield.py`: exact arithmetic in Q(sqrt 2, sqrt 3).
- `python/boundary/build_records.py`: builds `data/boundary/records.json` from the sources and
  the exact constructions. Deterministic; a test rebuilds and compares.
- `data/boundary/sources/research_round2.json`: the research prototype's round-2 output,
  imported with its numbers unchanged. `data/boundary/sources/search_2026-10-25.json`: the
  search's 11 configurations with their calls (method, seed, evaluations, date, value).
- `python/boundary/test_boundary.py`: tests; `python/boundary/fixtures/`: adversarial cases and
  the engine cross-check.
- `packages/boundary/scripts/check-records.ts`, `packages/boundary/scripts/gen-crosscheck.ts`,
  `packages/boundary/scripts/search-sources.ts`.
- The other containers (section 4A): `python/boundary/general.py` (the general exact certifier), `python/boundary/test_general.py`
  (its tests), `python/boundary/build_general_packings.py` (certifies the test packings), `packages/boundary/src/containers.ts` (the vertex lists),
  `packages/boundary/src/bounds.ts` (facts, lower bounds, check-size estimates), `packages/boundary/src/domain.ts` (the 152 problems),
  `packages/boundary/scripts/gen-containers.ts`, `gen-general-crosscheck.ts` and `gen-general-packings.ts` (the fixtures
  `containers.json`, `general_crosscheck.json`, `general_packings_raw.json` of `python/boundary/fixtures`).
- The Platonic solids (section 4C): `packages/boundary/src/pieces.ts` (the solids), `validate.ts` and `parts.ts` (the validator), `bounds.ts` (`pieceBounds`, `BOUNDARY_SOLID_MS`),
  `search/solids.ts`, `search/rot.ts` and `kernel-solids.ts` (the model, the rotations, the kernel), `python/boundary/platonic.py` and `qmulti.py` (the exact certifier and its field),
  `python/boundary/test_platonic.py`, `packages/boundary/test/solids.test.ts` and `solids-search.test.ts` (tests), `python/boundary/build_platonic_packings.py` (certifies the test packings),
  `packages/boundary/scripts/gen-platonic-packings.ts` and `gen-platonic-crosscheck.ts` (the fixtures `platonic_packings_raw.json` and `platonic_crosscheck.json`).
- The d-cubes (section 4B): `packages/boundary/src/ncube/` (`lattice.ts` words and the reach search, `validate.ts`, `bounds.ts`, `move.ts`, `model.ts`, `kernel.ts`, `seeds.ts`, `rule.ts`, `search.ts`),
  `python/boundary/ncube.py` (the exact certifier), `python/boundary/test_ncube.py` and `packages/boundary/test/ncube.test.ts` (tests), `python/boundary/build_ncube_packings.py` (the test packings),
  `packages/boundary/scripts/gen-ncube-crosscheck.ts` and `gen-ncube-bounds.ts` (the fixtures `ncube_crosscheck.json`, `ncube_bounds.json`).

```
python3 python/boundary/certify.py data/boundary/records.json
python3 python/boundary/test_boundary.py
python3 python/boundary/test_general.py
python3 python/boundary/test_ncube.py
python3 python/boundary/test_platonic.py
npx tsx packages/boundary/scripts/search-sources.ts
python3 python/boundary/build_records.py
npx tsx packages/boundary/scripts/check-records.ts
npx tsx packages/boundary/scripts/gen-crosscheck.ts
npx tsx packages/boundary/scripts/gen-containers.ts
npx tsx packages/boundary/scripts/gen-general-crosscheck.ts
npx tsx packages/boundary/scripts/gen-general-packings.ts
python3 python/boundary/build_general_packings.py
python3 python/boundary/build_ncube_packings.py
npx tsx packages/boundary/scripts/gen-platonic-packings.ts
python3 python/boundary/build_platonic_packings.py
npx tsx packages/boundary/scripts/gen-platonic-crosscheck.ts
npx tsx packages/boundary/scripts/gen-ncube-crosscheck.ts
npx tsx packages/boundary/scripts/gen-ncube-bounds.ts
```

Timing (Apple Silicon, shared machine): the 11 searches about 2.5 minutes in all (the longest,
balls P = 13 and 15, 44 and 49 s); building the records file 3 s; certifying all 32 records with
their exact constructions 2.3 s; the Python tests about 26 s (`test_general.py` 90 s, `test_platonic.py` 100 s); the engine check 0.3 s for 64
validations.

## 11. Sources

- E. Miller, I. Pak, Metric combinatorics of convex polyhedra: cut loci and nonoverlapping
  unfoldings, Discrete Comput. Geom. 39 (2008) 339-388, https://arxiv.org/abs/math/0312253
  (full text, sections 1-3, read by the research note): boundaries of convex polytopes in every
  dimension; flat and warped points; Prop. 1.2 (shortest paths avoid warped points); the
  exponential map. No packing.
- S. Devadoss, M. Harvey, Unfoldings and nets of regular polytopes (2021),
  https://arxiv.org/abs/2111.01359, and S. Devadoss, M. Harvey, D. Richter, Unfolding overlaps of
  the exceptional regular polytopes (2026), https://arxiv.org/abs/2608.30285 (abstracts and
  introductions, as read by the research note): every ridge unfolding of the tesseract is a net;
  nets, not packings.
- F. Buekenhout, M. Parker, Discrete Math. 186 (1998) 69-94: the count of the tesseract's nets,
  261 (as cited by Devadoss and Harvey; not opened).
- polyhedrainpolyhedra.wiki (A. Zarzuelo Urdiales, github.com/alejandrozu/platonic-packing): the
  cube-in-cube table used for the cube nofold value at two per cell, as cited by the research
  note (site and repository audited for PolyPack on 2026-10-23).
- PolyPack research note "3D pieces on the boundary of a 4-polytope", 2026-10-24: definitions,
  the lattice-word developments, the round-2 searches and their configurations.
- PolyPack boundary search (`packages/boundary/src/search`), runs of 2026-10-25: the 11
  improved configurations.
- [SURFACE.md](SURFACE.md) section 2.3: the distance form for disks, which section 4.2 lifts one
  dimension up.

## 12. On the site: the API and its rule

The server lists the space as `boundary` (title "Polytope boundaries", viewer
`/spaces/boundary.html`) through the same routes as the other spaces; the API page, `/api`,
describes them all. Everything below goes through the spaces registry, `packages/spaces`.

- `GET /api/spaces`: the space with its conventions and its problems (its `domain`: polytope, piece, P; the pieces and the P range depend on the polytope): 816 on the 4-polytopes (4C.1) and the d-cubes' 4,288.
- `GET /api/spaces/boundary/records`: every problem with a research record (the tesseract's 32) with `fields` `{polytope, piece, P}`, the
  research record (`certified`: score_cert, display, status), the lower bound (`bound`), the site's
  best entry and `best_score`. Boundary problems also carry `extra`:
  `{nofold: {value, label}, fold_gain}`, the nofold value of section 6 with its source sentence, and
  the record's fold_gain = 1 - score_float / nofold (section 7). Other spaces send no `extra`.
- `GET /api/spaces/boundary/problems/boundary/tesseract/ball/9`: one problem, with
  `certified_config` (the record's configuration at s = score_cert) and `site_config`. A problem of
  the other containers, `.../problems/boundary/24cell/ball/5`, answers with `status: "no_record"`, its lower bound (4A.4), no
  `certified` and no `certified_config`.
- `GET /api/spaces/targets?spaces=boundary`, until it was removed with the volunteers' work list on
  2026-10-29: the open problems with room above their bound, 28 of the 32. Left out: P = 1, and the
  records that meet their bound (balls P = 2 at sqrt 2, cubes P = 8 at 1). The volunteer client now
  draws among all 32 at random ([VOLUNTEER.md](../VOLUNTEER.md#the-end-of-the-work-list)).
- `POST /api/spaces/submit`: `{"space": "boundary", "id": ..., "config": {s, pieces}, "name": ...}`.
- `GET /api/spaces/entries/<id>`: a kept entry, with its configuration.

**The server's rule** for a configuration:

1. The id must be a problem of the domain, in its one written form (no leading zeros): balls, aligned cubes and the four other solids (rotated) on the tesseract, P = 1 to 16, balls and the five solids on
   the 5-cell, 16-cell, 24-cell, 120-cell and 600-cell, P = 1 to 24. P is at most 64 for the space and at most the
   catalogue's largest P (24) for a submission.
2. The configuration must be well formed, or the answer is 400 with what is wrong: s a finite
   number in (0, 1e4]; `pieces` a list; each piece `{cell, u}` with `cell` an integer from 0 to C - 1 (C = 8 on the tesseract, 5, 16, 24, 120, 600 on the others)
   and `u` three finite numbers, each at most 1e4 in absolute value; a cube or another solid may have `q`, four finite
   numbers whose norm is within 1e-9 of 1 (absent: no rotation); a ball has no `q`. Unknown fields are dropped.
3. A wrong number of pieces, or an s below the problem's lower bound (less a relative 1e-9), is an
   invalid packing: 422 with the reason, before any geometry.
4. `validateBoundary` (section 4.5's engine validator) at the touching tolerance 1e-9, with at most
   1024 chains of cells per piece: 422 with up to 20 reasons when it refuses.
5. An entry below score_cert by more than a relative 1e-9 must also pass at tolerance 1e-11
   ("beats the research record only within the float tolerance" otherwise). Every record passes it
   at its stored s. On a problem without a record (all but the tesseract's 32) every entry must pass it
   (`SPACE_STRICT_NO_RECORD_REASON`).
6. It is kept only when it beats the problem's best (record or site leader) by more than a relative
   1e-9, and it is stored as sent, with its known fields.

A site entry is checked in floating point and is never called certified. The research records are
the only certified values (section 4).

[CALC] Work: every research record, and every one of 1,449 random valid single balls and cubes,
needs at most 79 chains per piece, so the cap of 1024 refuses no real packing. Checks of 64 pieces
took 0.1 to 20 ms on an Apple Silicon laptop: the validator's volume bound and in-cell edge tests
refuse most bad input before any chain is enumerated. The server keeps the 32 records in about
55 kB of heap. On the other containers the test packings need at most 399 chains per piece (the 120-cell; 245 on the 5-cell and 600-cell),
and a check of 24 balls took 4 to 18 ms (5-cell, 16-cell, 24-cell, 600-cell) and 91 ms (120-cell) on loose layouts; balls on the 120-cell from P = 5
and on the 24-cell from P = 22 wait in the server's heavy-check queue (4A.6).

[CALC] Platonic solids: a check costs 1.2 to 4.4 ms a piece when every solid touches its neighbours, at most 2.2 ms a piece alone, and a solid reaches at most 813 chains of cells (400 random poses of each pair of container and solid, at 1 to 3 times the one-piece bound: the tetrahedron on the 120-cell where nothing fits; 47 to 566 for the others) and about 80 per piece in the test packings (4C.6), so the cap refuses no real packing; the server prices them at `BOUNDARY_SOLID_MS` per piece (4C.2) and queues those of 50 ms or more. The work caps: 1024 chains per piece, 400,000 tests of pairs of parts.

**Searches.** The board and the volunteer client run the same job, `runSpaceJob` in
`packages/spaces/src/search`, under the method name "boundary compression search": the boundary
package's `searchBoundary`, seeded, with an evaluation budget, from the best so far. The volunteer
client searches the space only when asked to (`--dimensions boundary`, or `all`); see
[VOLUNTEER.md](../VOLUNTEER.md).

```
curl -s --compressed https://polypack.fly.dev/api/spaces/boundary/records
curl -s --compressed https://polypack.fly.dev/api/spaces/boundary/problems/boundary/tesseract/cube/12
```

## 13. On the website: the board and the viewer

**On the board.** Choose **3D on 4D** under Dimension (its one kind, Polytope boundaries, is chosen
with it), or follow a link such as `/#space=boundary/tesseract/ball/9`; `/#space=boundary` opens
that problem. The Problem section offers Polytope (Tesseract, 5-cell, 16-cell, 24-cell, 120-cell, 600-cell), Piece (Balls, Cubes, Tetrahedra, Octahedra, Dodecahedra, Icosahedra) and Pieces (P),
each with the values the catalogue has. The stage shows the drawing (below), the panel the usual
values (the edge s of the packing shown, the certified value, the status, the lower bound, the
site best, whether the problem takes entries) and two of this kind:

- **Without folding**: the nofold value of section 6, with its derivation under "Why".
- **Folding gains**: 1 - best / nofold, where best is the problem's best edge (the research record
  or a better site entry), to one decimal of a percent. For a record this is section 7's gain, to
  within the 1e-12 of certification; a site entry that beats the record raises it.

Run optimizer runs section 12's search in the browser, one complete search call per round, of
300,000 energy evaluations (up to three times that while rounds get nowhere), about 1.5 to 2.2 s
each at P = 16 (balls and aligned cubes on the tesseract; on the other containers and for the Platonic solids a round is about two seconds of the model, at least 200 evaluations, and a solid's first round 2,000 + 1,000 P, which holds the search's first offer: 4C.6), and Submit sends a result that passes the server's rule
and beats the problem's best by more than a relative 1e-9; the board's notes
([BOARD.md](BOARD.md)) describe both. Pieces cannot be moved by hand.

**The drawing** is the three.js view of `packages/view3d`, the same in both places, drawn only
from the engine's parts (`pieceParts`), so what is seen is what the validator checks:

- **Cells**: the 8 cells side by side, each in its own coordinates; a face's colour names the cell
  glued to it there, so a piece that folds across a 2-face shows in two tiles in one colour.
- **Net**: the cells as a net (the Dali cross); **Fold** turns each cell about its 2-face into the
  fourth dimension, seen from a 4D eye, until at full fold it is the classic picture of the
  tesseract, a small cube inside a large one.
- **Slice**: the boundary cut by the hyperplane x_k = t, for an axis x, y, z or w and a position
  t; inside, the slice is the surface of a cube with each piece's cross-section on it, and at
  either end it is a whole cell. A strip of nine slices shows the sweep at once.

Drag to turn, or focus the drawing and use the arrow keys (Shift for larger steps), + and - to
zoom, Page Up and Page Down to follow one piece (its cells and centre are named below the
drawing), Escape to let it go, Home to reset. Choosing another P keeps the view, the fold and the
slice; another piece opens a new drawing in the view last used. three.js (166 kB gzipped) is
loaded only when a 3D on 4D problem is drawn. Without WebGL the drawing says so in a sentence; the
values, the optimizer and Submit still work.

**The viewer**, `/spaces/boundary.html`, is a page of the site built from
`apps/web/spaces/boundary.html` and `apps/web/src/boundary-viewer.ts`. It reads the catalogue
(`GET /api/spaces/boundary/records`) and opens any problem with the same picker, plus a Packing
menu that shows either the research record or the site's best entry. It shows the drawing, the same
values, and links to the problem on the board, to these notes, to the records table and to the
problem as JSON. `?id=boundary/tesseract/cube/12` or `#id=boundary/tesseract/cube/12` opens a
problem (the hash wins); the address follows the problem shown, and a problem that is not listed
opens the nearest P, with a sentence saying so.

**The dimensions page** files the space under 3D on 4D, with Records, Board, Viewer and Notes
links. Its records table starts with the balls and adds the columns Without folding and Folding
gains.

## 14. Kernel for the generic methods

`boundaryKernel({ polytope, piece, P })` (`packages/boundary/src/kernel.ts`) lets the generic
search methods ([docs/methods/dimensions/README.md](../methods/dimensions/README.md)) run here.
It reuses the boundary search's own parts, so a method sees the same model as section 12's search.

- **State.** Each piece has an anchor (a cell and a centre, kept in the state's `aux`) and a
  displacement x from it, 3 lengths in the anchor's frame. The piece's centre is where the
  straight path from the anchor along x ends, unfolded across every 2-face it crosses (the
  search's charts). `rebase` moves the anchors to the centres. The size is s, the score. A path
  that meets an edge, or one longer than 64 edges, leaves its piece at the anchor. That length is
  far beyond the boundary's diameter [CALC: over 80,000 random pairs of points the shortest
  geodesic was at most 2.3 s].
- **Moves.** A move changes only the piece's x, so writing x back undoes it. `jiggle` adds to x.
  `randomPose` draws a uniform target (a cell and a point in it) and writes the displacement from
  the anchor to it along the shortest geodesic the lattice words give (section 4.1).
- **Energy.** The search's penalty energy: lattice words for balls, aligned boxes for cubes
  (sections 4.2 and 4.3 state the exact rules they follow). It takes per-piece weights, can be
  restricted to one piece's terms, and has an exact gradient in x and in s [CALC: checked against
  central differences in `packages/boundary/test/kernel.test.ts`].
- **finish** is the search's last stage with one addition: the scale is bisected down as well as
  up, so a state with room is drawn in until its pieces touch. The result is then nudged up until
  `validateBoundary` accepts it at tolerance 0.
- **validate** is the server's rule of section 12, rebuilt in `kernel-validate.ts` because this
  package does not depend on the spaces registry. It checks P (at most 64), the form in the
  server's words, s against the problem's lower bound less 1e-9 relative, and then runs
  `validateBoundary` at tolerance 1e-9 with at most 1024 chains per piece. The bound is the one
  the records store: for each of the two bounds of section 5, the largest binary64 value not
  above its exact value, decided in integer arithmetic. The score is s whenever the
  configuration is well formed, as on the server. Only the strict check at 1e-11 of an entry that
  beats a research record is left out, since it needs the record. Configurations from `finish`
  pass it anyway. `packages/spaces/test/boundary-kernel.test.ts` compares the verdict, the score
  and the reasons with the server's `validateSpaceSubmission` on every record, on 26
  variants of each, and on finished configurations.
- **fromConfig** takes only what the server could take: a well-formed configuration within its
  limits (s and coordinates at most 1e4). A start outside them gives no state, so no method works
  on an absurd one. For the same reason there is no kernel for P above 64.
- **Cubes stay axis-aligned.** The cube model's exact rule covers aligned cubes only, and the
  kernel produces nothing its model cannot judge. A start with a turned cube is refused, unless
  the turn is a symmetry of the cube. For the same reason the kernel's `bound` for cubes, where
  a method may stop, is max(1, (P/8)^(1/3)): no aligned cube fits below s = 1. `validate` uses the
  problem's own bound, which allows turned cubes. Balls use the proven bound of section 5 for both.
- **Random starts.** They are the search's random layouts. Above P = 16, outside the catalogue,
  cubes start at an edge where they fill at most a quarter of the boundary, each inside its cell
  [CALC: at edge 2, random starts of P = 26, 27 and 31 to 64 cubes often relaxed to nothing a
  finish could use].
- **Capabilities.** `seeds`: the nofold layouts of section 6 (n pieces per cell at the centre, on
  a diagonal or on a grid) and, for cubes, the search's two rings. `shape` for balls only, where
  it changes nothing (a ball is its own inscribed ball). Cubes have no soft shapes, because the
  cube model has no rounded cubes. No piece has an orientation, so there is no `swap`.
- **Costs** [CALC, in process CPU time on an Apple Silicon laptop]. One energy() call takes 2 to
  15 microseconds at P = 4 to 16. Balls cost most below s = 2, where the longer word table is
  needed. A finish is one validator call at tolerance 0 plus 20 to 45 evaluations of the model's
  rule; a validate is one validator call. On squeezed random starts and on the records at
  P = 4, 8, 12 and 16, a finish cost a median of 930 energy() calls and a validate 960, each
  ranging from about 360 to 4300. Cubes at P = 4, near s = 1, cost the most. The kernel states
  `finishCost = validateCost = 1000`.
- **Records.** Every record converts, passes the server's rule, and finishes about 1e-12 below its
  score_cert. That gap is the certification's own scaling (section 4.4), not an improvement.
- **The other containers** (balls only; 4A.6). The energy is the chain metric's (`ChainMetric`, chains of cells enumerated from each centre at
  the ball's diameter), with the same exact gradient in x and in s, the weights and the focus; `dE/ds` holds the centres fixed and scales the
  developments' translations and the nearest hinge points with s. `randomPose` develops a uniform target through the chains at a growing reach
  (at most `POSE_CHAINS` = 1200 chains), and, for a target too far for that (almost always on the 600-cell), uses the end of a random geodesic of
  up to 2.5 edges instead. The seeds are balls at the centres of P cells (valid from s0 of 4A.4), or a sequential random layout when P exceeds the
  number of cells. The bound is 4A.4's, the form check takes cells 0 to C - 1, `finish` and `validate` are unchanged. [CALC] One energy evaluation
  with its gradient costs 0.1 to 0.8 ms per ball (4A.6).
- **The Platonic solids** (4C.6; `packages/boundary/src/kernel-solids.ts`). `dim` is 6: the displacement (3 lengths) and a rotation vector (3 radians) about the anchor's orientation. The orientation at the piece's centre is M exp(w) R_a,
  M the linear part of the gluings crossed (a piece whose path meets an edge keeps its rotation), and the gradient in w is Jl(w)^T det(M) M^T times the model's gradient in the rotation increment. A typical step is the solid's inradius and
  inradius / circumradius radians. `swap` gives each of two pieces the other's orientation; every Platonic solid is mirror-symmetric, so an orientation of the wrong handedness for the walk is replaced by its mirror image. `seeds`: the cells
  themselves when the solid is the container's cell and P is at most the cell count, else a sequential random layout accepted piece by piece by the model; no `shape`. `finishCost` is 60 and `validateCost` 15 energy evaluations
  [CALC: a finish took 12 to 55 energy evaluations, a validator call 0.1 to 1 on random states and up to about 15 on a crowd of touching solids]. `bound` is 4C.4's. The tesseract's `cube` keeps the aligned-cube kernel above.
