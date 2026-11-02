# 3D on 4D, polytope boundaries: balls (diameter 1) and cubes on the boundary of a regular 4-polytope

**Dimension.** Under Dimensions this is **3D on 4D**: solid 3D pieces packed on the 3-dimensional
boundary of a 4D polytope, with no container around them. Problem ids and records use the space
name `boundary`.

P balls of diameter 1 (radius 1/2) or P unit cubes (edge 1) are packed on the boundary of a tesseract of
edge s. Since 2026-11-01 balls can also be packed on the boundary of the other five regular
4-polytopes, the 5-cell, 16-cell, 24-cell, 120-cell and 600-cell (section 4A). Also since 2026-11-01 balls and
axis-aligned cubes (of dimension d - 1) can be packed on the boundary of the d-cube, d = 5 to 10 (4D on 5D to 9D on 10D, section 4B), and, since 2026-11-02, balls on the boundary of the d-simplex, d-orthoplex and d-demicube for the same d (section 4C). The score is the smallest s,
the polytope's edge. Since 2026-10-27 balls have diameter 1; before, radius 1, and every ball score was twice as large (cube scores are unchanged). The sections below are stated in the site's units; the certifier and the ball code may work in radius-1 units internally (a factor of 2, exact) and convert at their edges. Every record is certified by an exact checker
(`python/boundary/certify.py`, exact arithmetic only), so its validity is proven, not just
tested. Records: `data/boundary/records.json`, from the PolyPack research prototype (2026-10-24),
the PolyPack boundary search (2026-10-25) and exact constructions (section 4.6). They are the tesseract's. The five other containers have
problems and lower bounds but no records (section 4A.7); their test packings are certified by the
general exact certifier (`python/boundary/general.py`, section 4A).

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

There are 32 problems: P = 1..16 for balls and for cubes.

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
exact certifier is `python/boundary/general.py`. Only balls are covered here (4A.8 says what the Platonic pieces need).

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
  its 32 problems, P = 1 to 16, ball and cube). A cube on any of the five is not a problem. 152 problems in all.

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

- **Problems.** 152 in all: the tesseract's 32 (P = 1 to 16, ball and cube), and balls on the 5-cell, 16-cell, 24-cell, 120-cell and 600-cell, P = 1 to 24 each
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

- **Platonic pieces are not built** (the tetrahedron, cube, octahedron, dodecahedron and icosahedron, rotated, on the six containers). What
  remains, in the order it would be done: (1) the float validator's parts for these shapes (`parts.ts` knows the unit cube only; the clipping
  and the separating-axis code of `polyhedron.ts` are general, and `validateCubes` would become a validator of polytope pieces, with the hinge
  depth inside a piece taken over its faces instead of a cube's axes); (2) a penalty model for rotated polytope pieces in the search and the
  kernel (today's cube model, `search/cubes.ts`, is axis-aligned: the state needs a rotation per piece carried across ridges by the gluings'
  linear parts, the energy is the separating-axis depth of the pieces' parts in each cell, and the kernel's `dim` becomes 6 with a `swap`
  for equal pieces); (3) the exact certifier's polytope-part test: the piece, turned by an exact rational rotation (a rational quaternion gives
  one), intersected with each developed cell and window, and the parts of different pieces, and of one piece, compared in each cell by an exact
  separating-axis test, with the hinges outside the piece's interior; (4) the field: in base units a piece has edge lambda = e0 / s times the
  unit-edge solid's coordinates, which involve sqrt 2 (tetrahedron, octahedron) or sqrt 5 (dodecahedron, icosahedron); with e0 = sqrt 2
  (5-cell, 16-cell, 24-cell) the tetrahedron and octahedron come out rational and the other three in Q(sqrt 2) or Q(sqrt 2, sqrt 5), and
  on the 120-cell and 600-cell the tetrahedron and octahedron need Q(sqrt 2, sqrt 5), a field of degree 4; (5) test packings of each piece on
  each container: the cells themselves (P = C pieces at s = 1 for the tetrahedra of the 5-cell and the dodecahedra of the 120-cell, the
  tesseract's cubes) are exact ones. Nothing of this is in the code.
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
- Not covered here: the simplex, orthoplex and demicube boundaries for d >= 5, whose dihedral angles are not right angles. Balls on them are section 4C (2026-11-02); polytope pieces on them are not built.


## 4C. Balls on the d-simplex, d-orthoplex and d-demicube: (d-1)D on dD, d = 5 to 10

Added 2026-11-02. Under Dimensions each of **4D on 5D**, **5D on 6D**, ..., **9D on 10D** now holds four containers: the d-cube of section 4B and the
**d-simplex**, the **d-orthoplex** and the **d-demicube**, whose boundaries are made of (d-1)-simplices, or of (d-1)-simplices and (d-1)-demicubes,
meeting at angles that are not right angles. Pieces are balls of diameter 1 only; polytope pieces on these boundaries are a later phase (4C.10).
The rule is 4.2's, which is dimension-free; what changes is how the shortest paths are found. Section 4B's lattice is gone (only the cube has one), so the engine is
the chain development of 4A taken to R^(d-1), with one change that makes it cheap: **pairs and loops are found by joining two half-chains of
reach 1/2** (4C.2). Code: `packages/boundary/src/nd/`; exact certifier: `python/boundary/ndpoly.py`. Nothing is shown in Tighter.

### 4C.1 The containers

Each container has unit edge scaled by the score s, and is built (`nd/container.ts`; `ndpoly.py` repeats it in exact arithmetic) from the classical
coordinates of edge sqrt 2 [EST: Coxeter, *Regular Polytopes*, 7.2 and 8.7 (the demicube); the facet counts below are the literature's], divided by sqrt 2:

| container (d = 5 ... 10) | vertices | cells (facets) | ridges | cell types | cos of the dihedral angle |
|---|---|---|---|---|---|
| d-simplex | e_0 ... e_d in R^(d+1) | d + 1 = 6, 7, 8, 9, 10, 11 | C(d+1, 2) = 15, 21, 28, 36, 45, 55 | (d-1)-simplex | 1/d |
| d-orthoplex | +-e_i in R^d | 2^d = 32, 64, 128, 256, 512, 1024 | d 2^(d-1) = 80, 192, 448, 1024, 2304, 5120 | (d-1)-simplex | (d-2)/d |
| d-demicube | (+-1/2)^d, an even number of minus signs | 2d + 2^(d-1) = 26, 44, 78, 144, 274, 532 | 2^(d-1) d + 2d (d-1) = 120, 252, 532, 1136, 2448, 5300 | 2d (d-1)-demicubes and 2^(d-1) (d-1)-simplices | 0 (demicube, demicube); 1/sqrt d (demicube, simplex) |

The simplex lies in the hyperplane sum x = 1 of R^(d+1) (its coordinates in the engine are the barycentric ones), the others in R^d. The 5-demicube's facets are 16 5-cells (4-simplices)
and 10 16-cells (4-demicubes); in general the demicube's (d-1)-simplex facets are all surrounded by demicube facets (a simplex facet has d ridges, each shared with a
demicube facet), and two simplex facets never share a ridge. The simplex-to-demicube dihedral cosine 1/sqrt d is irrational for d != 9 (4C.2, 4C.4). [CALC; the tests read the vertex set by brute force and compare the cells, their neighbours and their hinges for d = 5, 6, 7.]

- **Cells.** Simplex: cell i is the facet x_i = 0 (vertices e_j, j != i); normal (1 - (d+1) e_i) / sqrt(d (d+1)), offset 1 / sqrt(d (d+1)) (base units). Orthoplex: cell c is the facet
  sigma . x = 1 with sigma_i = -1 exactly when bit i of c is set (vertices sigma_i e_i); cells c and c xor 2^k are adjacent. Demicube: cell 2 i + (a > 0) is the facet x_i = a/2 (a = +-1),
  adjacent to the other demicube facets x_j = b/2 (j != i) and to the simplex facets with w_i = a; cell 2d + k is the simplex facet w . x = (d-2)/2 at the k-th sign vector w of odd parity (in increasing order of its
  bit mask, bit i set exactly when w_i = -1), adjacent to its d demicube facets x_i = w_i/2.
- **Hinges.** The singular set is the (d-3)-faces. A hinge of a cell is the intersection of two of its ridges that has dimension d - 3 (`hinges` in `container.ts`): every pair of ridges of a simplex or
  orthoplex cell; of a demicube cell, a pair of demicube-facet ridges x_j = b, x_k = b' (j != k, d >= 6) or a demicube-facet ridge x_j = b and a simplex-facet ridge with w_j = b; of a simplex cell, every pair.
- **Local coordinates.** u has n = d - 1 numbers: x = s o_c + sum_k u_k f_k, o_c the centroid of the cell's vertices and f_0, ..., f_(n-1) the cell's frame: Gram-Schmidt (twice, for accuracy) on the projections of e_0, e_1, ... onto the cell's
  tangent space, skipping a vector whose residual is under 1e-9, until n are found. On a demicube facet x_i = a/2 this is the other coordinates in increasing order. The engine's ambient coordinates are R^D (D = d + 1 for the simplex, d otherwise).
- **Ids.** `boundary/<d><simplex|orthoplex|demicube>/ball/<P>`, d = 5 ... 10, P = 1 to 24, no leading zeros: 18 x 24 = 432 problems (`boundary/5simplex/ball/1` ... `boundary/10demicube/ball/24`). Cubes are not a problem on these containers.
  Python's `ID_RE` (`certify.py`) names them and keeps its `\Z`.

### 4C.2 Developments, windows and joins

**The development map** [CALC]. Let c and c' be adjacent cells with outward unit normals n and n' and offsets h, h' (n . x <= h on the polytope), and q a point of their ridge. The reflection
    sigma(x) = x - 2 (m . (x - q)) / (m . m) m,   m = n + n',
is an involution and an isometry of R^D; it fixes the ridge pointwise (m . (x - q) = 0 there); it maps the hyperplane H_c' onto H_c; and it lays c' on the *far side* of the ridge, because for x on c'
    n' . sigma(x) = h + h' - n . x >= h',    (with equality exactly on the ridge)
(use n . m = n' . m = (m . m)/2, true for unit normals), i.e. sigma(c') is outside the half-space of c's facet c'. (The other reflection through the ridge, with normal n - n', folds c' back over c.) On H_c, sigma coincides
with the rotation of the plane (n, n') that takes n to n' (4A.3's gluing); the tests compare them. Because sigma is an affine reflection, an exact development is a product of such maps, with no rotation matrix of d - 1 dimensions to write:
**with the containers' coordinates every reflection is rational for the simplex, the orthoplex and the demicube's ridges between two demicube facets, and for the demicube's ridges between a demicube facet and a simplex facet exactly when d = 9**;
for d = 5, 6, 7, 8, 10 the cosine 1/sqrt d of that dihedral angle is irrational, sigma(x) of a rational x is not rational, and the exact certifier works in Q(sqrt k), k the squarefree part of d (4C.4). [CALC, checked exactly on all d = 5 ... 10 (`test_ndpoly.py Developments.test_rationality`); the first guess, that every reflection is rational, is true for two of the three families.]

**Half-chains.** From a centre u of cell c_0, a node is a cell c reached by a sequence of ridges, with a = T^(-1) u, the source as seen in the hyperplane of c (the development T of c into the hyperplane of c_0 is the composition of the reflections; a child's a
is sigma(a), sigma an involution). A ray from u crosses the ridge of the neighbour g of c only if u is strictly inside the half-space of g as seen from the developed cell: S_g = h_g - n_g . a > 0 (all in the scale s), and then the path is at least as long as the distance from a to the ridge, which is at least
    sqrt( S_g^2 / nu_g^2 + max_k (viol_k / tn_k)^2 ),
where nu_g^2 = 1 - (n_g . n_c)^2, the foot x_p of the perpendicular from a on the ridge's flat has viol_k = n_k . x_p - s h_k > 0 for the cell's other facet k, and tn_k is the length of that facet's normal tangent to the flat. [CALC] The first term is the distance from a to the ridge's flat inside H_c; inside that flat, the nearest point of the ridge lies at least viol_k / tn_k from x_p
(any single facet of the cell is a relaxation of the ridge), and the two are orthogonal. A node is kept when this bound is at most the reach; so a kept-node set contains every chain a straight path of length at most the reach can follow. **No window is stored or clipped**; the windows are tested where they are needed.

**Windows.** For a step (cell p, crossed neighbour f, S_k = s h_k - n_k . a_p) a ray from u with direction z (in p's frame) crosses the plane of f at tau = S_f / (n_f . z), needs n_f . z >= 0, and at the crossing satisfies the facet k of p: S_f (n_k . z) - S_k (n_f . z) <= 0 for every other neighbour k of p
(multiply n_k . (a_p + tau z) <= s h_k by n_f . z). The window of a chain is the cone of the directions satisfying these for every step, carrying z back from the end cell's frame by the reflections' linear parts (z_j = R_j z_(j+1)). **Lemma 1** [CALC, as 4A.2]: a ray in the closed window crosses the ridges of the chain in order, each inside the closed ridge,
and runs inside each developed cell between consecutive crossings (each developed cell lies on the far side of the plane of the ridge it was entered through, a ray crosses a plane once).

**Joins.** *Lemma 2* [CALC]. Let C1 and C2 be chains from u1 and u2 ending in one cell c, a and b the images of u1 and u2 in the hyperplane of c, and w = b - a. If w is in the closed window of C1 and beyond its last ridge (n_f . z >= S_f for its last step), and -w in the closed window of C2 and beyond its last ridge,
then [a, b] is the development of a path of length |w| from u1 to u2. Proof: by Lemma 1 the ray from a through w crosses the ridges of C1 in order inside the closed ridges, its crossing p of the last one at a parameter <= 1; likewise from b, the crossing p'. Both lie on the boundary of the convex cell c, so [p, p'] lies in c; and p comes before p' on [a, b], for p' is in c, which lies in the closed half-space of C1's last ridge, while [a, p) lies strictly outside it. So a, p, p', b are in this order and the segments
map back to a path of the boundary. *Completeness:* a shortest path of length L < 1 from u1 to u2 avoids the hinges [EST: Miller, Pak 2008, Prop. 1.2], so it is a straight segment of a development; its midpoint m is at distance L/2 < 1/2 from both ends, in a cell c; the halves are straight hinge-avoiding paths of length < 1/2 from u1 and from u2, so by 4A.2's argument (applied to the half-length) the chains
they follow are kept, and the join of those two chains is found with |w| = L. *Hence pairs and loops need only the half-chains of reach 1/2 of each centre*: a loop is a join of two different half-chains of the same centre (a loop of length L is a path from u to u, with its midpoint at L/2 from both ends; the trivial pair, one chain with itself, has w = 0 and is skipped; a centre on a ridge has the chains of both its cells).
[CALC] This is the step that makes the work bearable: a ball at the centre of a facet of the 10-simplex (edge 5.8) has 12 half-chains of reach 1/2, while chains of reach 1 from a facet's centre number in the thousands at d = 7 already ([CALC, a prototype with exact window emptiness tests and exact gate distances, 2026-11-02]: chains of reach 1 from random points of a 7-simplex facet of edge 3.2 numbered 3,000 to 5,000 even after the gate-distance pruning, and 27,000 to 53,000 with only the windows; at reach 1/2 the same points had 8 to 31 (d = 7), 9 to 59 (d = 8), 11 to 100 (d = 10) chains.)

**Hinge distance.** The distance from u to the hinges is the least, over the nodes and the hinges of their cells, of the distance from a to the part of the hinge inside the node's window and cell: a convex quadratic program. The nearest point of {x : A x = b, C x <= e} to a is found by projecting on the flat of the equalities and, if the projection violates some inequality, trying each violated one as an
added equality (the nearest point lies on a violated constraint's plane, because the normal-cone vector of the polytope at the nearest point has a positive inner product with some active constraint's normal, which the projection violates), recursively; a branch is cut as soon as the distance of its projection plus the largest violation measured inside its flat exceeds the cap (the reach), and rows whose half-space contains the whole ball of the cap about a are dropped.
A hinge is looked at only when both its ridge planes and the flat they span (a 2x2 system) are within the reach. [CALC; checked against random walks: walking just past the nearest point meets a hinge, just short of it does not, d = 5, 6, 9.]

**Finiteness** [SPEC]. That the enumeration stops for every input is not proved. The reach bound, the cell's structure and the cap keep it small in practice (4C.8); a failure of finiteness could cost time, never soundness: a check that exceeds its step budget is refused.

### 4C.3 Rules

4.2 unchanged, for the (d-1)-ball of diameter 1: valid exactly when every centre lies in its closed anchor cell, is at geodesic distance >= 1/2 from the hinges, every two centres are >= 1 apart and every geodesic loop at a centre is >= 1 long. (The engine tests loops, as joins of two different half-chains of one centre; [CALC] none shorter than 1 occurs in any of the 414 cross-check cases or the 72 packings.) A centre within the snap distance eta = max(2^-52 max(1, s), min(1e-14 max(1, s), tol/10)) of a ridge plane is moved onto it (the distance moved is charged to the verdict) and belongs to both cells, as 4A.

### 4C.4 The exact certifier

`python/boundary/ndpoly.py`, reached through `certify.py` (`boundary/<d><family>/ball/<P>` ids). Exact only: Fractions for the simplex and orthoplex; for the demicube the field Q(sqrt k), k the squarefree part of d (5, 6, 7, 2 for 8, 1 for 9, 10), implemented as integers (n + m sqrt k)/q, sign decided exactly. It shares no code with the engine and repeats only its conventions (cell order, frames). The base polytope has edge sqrt 2 and outward normals N_c of one common squared length kappa (d(d+1), d, d), so no square root of a
norm is taken; a ball of diameter 1 on the polytope of edge s has diameter sqrt 2 / s in base units, and every test compares a squared base length with 2/s^2 times a rational.

- **Reading a configuration.** As 4A.3: the ambient point s o_c + sum u_k f_k (floats, as the engine computes it) is written as an affine combination of n + 1 independent vertices of the cell with dyadic weights (rounded to 2^-60, summing to 1), and the exact point is the same combination of the exact vertices; the cell's frame is recomputed in Python floats from the same Gram-Schmidt. A point outside its closed cell by less than 1e-12 s is moved onto the facet planes it violates; more is an error.
- **Reading discrepancy** [CALC]. The exact point differs from the engine's float point by at most **1.6e-16 s** (6 random centres per container, 108 in all, `Reading.test_the_exact_point_is_within_1e_15_of_the_engines`). Over the cross-check (`fixtures/nd_crosscheck.json`, 144 random layouts of 1 to 3 balls on all 18 containers) the engine's clearances, pair distances and loops differ from the exact ones by at most **9.4e-16** (386 values), and
  the same distances are found or not found (a distance is reported only below its threshold: 261 clearances below 1/2, 125 pair distances below 1, no loop below 1); the chain code accumulates a few more roundings than the reading, which is why the strict slack below is generous.
- **Strict validity** asks for room for those differences: a pair distance or a loop must exceed 1 by more than `READING_PAIR` = 2e-13 s, a clearance 1/2 by more than `READING_CLEAR` = 1e-13 s (250 to 600 times the measured discrepancy). The reach is widened to match, and the plain (non-strict) test uses reach exactly 1/2.
- **Fail closed.** A check that needs more than `MAX_STEPS` = 2,000,000 steps (nodes created, projections solved, join candidates) or more than `MAX_NODES` = 20,000 chains from one centre raises `NDError` ("the exact check needs more than ..."); `certify.py` turns it into a `CertError`, never a verdict.
- **Cost.** All 72 test packings are certified strictly in 31 s: at most 0.15 s each on the simplex and orthoplex, 0.7 to 3 s on the demicube for 24 balls, and 17 s for the slowest, one ball at the centre of an 8-demicube facet (5,391 steps; the big demicubes' hundreds of ridges per cell and Q(sqrt k) arithmetic make a step cost milliseconds).

### 4C.5 Lower bounds

Both bounds are proved, stored as the largest binary64 value not above the exact one, decided in integer arithmetic (`nd/bounds.ts`; the exact check is `ndpoly.exact_bounds`, which agrees with the engine on all 432 stored values, both directions: the value passes, the next float up fails).

- **Volume** [CALC]. The centre of a ball is >= 1/2 from the hinges, so the geodesic ball of radius 1/2 about it is flat, a Euclidean (d-1)-ball; the developments are injective on the interiors of the balls' parts and the balls have disjoint interiors, so
  P V_n 2^-n <= s^n vol_1, V_n the volume of the unit n-ball, n = d - 1, vol_1 the (d-1)-volume of the boundary at unit edge: a regular n-simplex has volume sqrt(n+1) / (n! 2^(n/2)), an n-demicube 2^(-n/2) (1 - 2^(n-1)/n!) (the n-cube of edge 1/sqrt 2 less 2^(n-1) corner simplices, each (1/sqrt 2)^n / n!);
  the simplex's boundary has d + 1 simplices, the orthoplex's 2^d, the demicube's 2^(d-1) simplices and 2d demicubes. pi and the square roots are rounded the safe way (pi down to 60 digits, sqrt d up).
- **One ball** [CALC, as 4A.4]. A centre needs geodesic distance >= 1/2 from the hinges, hence from those of its own cell; every point p of a cell lies within rho of a hinge of the cell, rho the largest distance from the cell's centre O to a hinge. Proof: p is on a segment from O to a point q of a facet F of the cell; q is within r_F of the boundary of F
  (the farthest point of a polytope with a symmetry group fixing its centre from its boundary is the centre; r_F its inradius), so within r_F of a facet e of F, a hinge of the cell; the distance to the convex set e is convex, so d(p, e) <= (1 - t) d(O, e) + t d(q, e) <= max(d(O, e), r_F) <= rho, as r_F <= d(O_F, e) <= d(O, e) (O_F the foot of O on F's plane). So s rho >= 1/2.
  The cells' hinges are all at s / sqrt(d (d-2)) from the centre of a simplex or orthoplex cell, and the (d-1)-simplices of the demicube at s / sqrt(n^2 - 1); a (d-1)-demicube's hinges between a demicube facet and a simplex facet at sqrt((n^2 - 5n + 8) / (8 (n-1))) s (n = d - 1; the demicube-demicube hinges, s/2, are nearer for n >= 6 but the bound uses the largest).
  So **s >= sqrt(d (d-2)) / 2** on the simplex and orthoplex: 1.936492, 2.449490, 2.958040, 3.464102, 3.968627, 4.472136 for d = 5 ... 10; and **s^2 >= 2 (n-1) / (n^2 - 5 n + 8)** on the demicube: 1.224745, 1.000000, 0.845154, 0.738549, 0.661438, 0.603023.
  *The demicube's bound is not tight for d >= 7.* The centre of a demicube facet has clearance exactly s/2 (the demicube-demicube hinges x_i = a, x_j = b are s/2 from it), so one ball fits from s = 1 for d >= 6 (and from sqrt 6 / 2 for d = 5, the same as the bound there), and a hill-climbing search (6,000 steps from random starts, d = 6 and 7) found no point of a facet with a larger
  clearance; whether the sharper bound s >= 1 holds for d >= 7 is open [SPEC]. The test packings of one ball therefore sit at s = 1 on the demicubes d >= 7 (the exact layout), and at the bound on the others.
- Nothing else is claimed. The 24-ball packings of 4C.7 are 1.7 to 2.5 times the bound; no optimality is claimed for any problem.

### 4C.6 Checks around the certifier

- **Engine against certifier.** `fixtures/nd_crosscheck.json` (`scripts/gen-nd-crosscheck.ts`): 144 seeded random layouts, 8 per container, of 1 to 3 balls (later balls placed at the end of a random geodesic of length 0.7 to 1.2 from an earlier one, so that most pairs cross ridges), with the engine's clearance, pair distances and loops; the exact ones agree to 9.4e-16 (386 values) and every verdict agrees.
- **Adversarial cases** (270, same file): per container, ball 2 at the end of a geodesic of length 1 + delta from ball 1 (delta = 1e-3, 1e-6, 1e-9, 1e-12, -1e-12, -1e-9, -1e-6), a ball moved away from its cell's centre until its clearance is 1/2 + delta (bisection on the engine's measure; delta = 1e-3, 1e-9, 1e-12, -1e-12, -1e-9, -1e-3), and a ball exactly on a ridge with a second one 0.9 and 1.05 away across it.
  The check: **the float validator at tolerance 0 never accepts what the exact certifier refuses** (an engine-accepted configuration has exact margin >= -2e-14, the reading error; none had a negative one), and an exact strict verdict implies the engine's acceptance. All 270 were run once by hand: the engine accepts 54, refuses 216, and the exact certifier finds 54 strictly valid, the same ones (161 s).
  The test (`Adversarial`) runs 214 of them, all but the 7- to 10-demicube pair and hinge cases, which are slow in exact arithmetic, in 14 s: 42 accepted, 172 refused, 42 strictly valid.
- **Development maps and walks** (`nd.test.ts`): the reflection is an isometry, an involution, fixes the ridge, lays the neighbour on the far side and agrees with the rotation; the join search finds, for 12 random geodesic walks per container (hinge-free, lengths 0.3 to 1), a path at most as long as the walk, and walking the join's direction ends at the other centre; walking the hinge search's nearest point just past it meets a hinge.
- **Known distances** (`nd.test.ts`): a cell's centre is s / sqrt(d(d-2)) from the hinges (simplex, orthoplex), s/2 (demicube facets, d >= 6) and s / sqrt 6 (5-demicube); two simplex cells' centres are s sqrt(2/(n(n+1))) apart; two adjacent demicube facets' s / sqrt 2.
- **Fail-closed**: a step budget cuts a check off with "work limit" for every budget below the steps it needs, and accepts at exactly that many; the exact certifier raises at its own limits.

### 4C.7 Test packings

No records: the 432 problems ship with lower bounds and test packings only (`python/boundary/fixtures/nd_packings.json`, 72 packings, 4 per container; the cheap constructions of this phase, no long run):

1. **P = 1**: one ball at the centre of a cell, at the edge from which it is valid: sqrt(d(d-2))/2 (simplex, orthoplex), sqrt 6/2 (5-demicube), 1 (the other demicubes), times 1 + 1e-9.
2. **The exact layout** [CALC: the derived edge is the smallest at which the layout is valid, certified exactly at 1 +- 1e-9 for every d = 5 ... 10]. Simplex: d + 1 balls at the centres of all the cells. Two cells are adjacent, and the centres are 2 r apart along the straight path across their ridge (the centres are mirror images in the ridge's hyperplane, each r from it), r = s / sqrt(2 n (n+1)) the inradius of the (d-1)-simplex of edge s: valid from
   **s = sqrt(d (d-1) / 2)** (3.162278, 3.872983, 4.582576, 5.291503, 6, 6.708204); the one-ball bound is below it. Orthoplex: two adjacent cells, from the same s. Demicube: 2d balls at the centres of the 2d demicube facets, valid from **s = sqrt 2**: two adjacent facets' centres are s/(2 sqrt 2) each from their right-angled ridge, s / sqrt 2 apart; opposite facets' centres are farther
   (through a third facet, sqrt 2 s). Below these edges the nearest pair overlaps (the exact certifier reports a path of length 1 - 1e-9 for 1 - 1e-9).
3. **P = 4 and P = 24**: short seeded runs of the boundary search (seed 7, 3,000 evaluations: the kernel's structured seeds, finished by its exact finishing, which draws the edge in until the balls touch; no optimiser), not good packings.

Each is scaled about the centres of the cells (4.4's argument: lam_k = 1 + k 1e-12) to the first k at which the exact certifier finds it **strictly valid with every margin at least 1e-12** (clearance over 1/2, pair distance and loop over 1; the search looks 2e-12 past the thresholds so that a smaller margin is seen and a larger one is certain): k = 0 in 36 packings, 1 in 12, 2 in 24. Shrunk by 1e-6 each is refused by the exact certifier and the engine.

| container | cells | bound, P = 1 | bound, P = 24 | P = 1: s | layout: P, s | P = 4: s | P = 24: s | s / bound at P = 24 |
|---|---|---|---|---|---|---|---|---|
| 5simplex | 6 | 1.936492 | 2.697733 | 1.936492 | 6, 3.162278 | 3.162278 | 5.000720 | 1.85 |
| 6simplex | 7 | 2.449490 | 2.746556 | 2.449490 | 7, 3.872983 | 3.872983 | 5.700797 | 2.08 |
| 7simplex | 8 | 2.958040 | 2.958040 | 2.958040 | 8, 4.582576 | 4.582576 | 6.813685 | 2.30 |
| 8simplex | 9 | 3.464102 | 3.464102 | 3.464102 | 9, 5.291503 | 5.291503 | 7.142530 | 2.06 |
| 9simplex | 10 | 3.968627 | 3.968627 | 3.968627 | 10, 6.000000 | 6.000000 | 8.280686 | 2.09 |
| 10simplex | 11 | 4.472136 | 4.472136 | 4.472136 | 11, 6.708204 | 6.708204 | 9.225848 | 2.06 |
| 5orthoplex | 32 | 1.936492 | 1.936492 | 1.936492 | 2, 3.162278 | 2.747523 | 3.567716 | 1.84 |
| 6orthoplex | 64 | 2.449490 | 2.449490 | 2.449490 | 2, 3.872983 | 3.595257 | 4.069306 | 1.66 |
| 7orthoplex | 128 | 2.958040 | 2.958040 | 2.958040 | 2, 4.582576 | 4.453174 | 4.975866 | 1.68 |
| 8orthoplex | 256 | 3.464102 | 3.464102 | 3.464102 | 2, 5.291503 | 5.715457 | 6.061795 | 1.75 |
| 9orthoplex | 512 | 3.968627 | 3.968627 | 3.968627 | 2, 6.000000 | 6.202058 | 7.327412 | 1.85 |
| 10orthoplex | 1024 | 4.472136 | 4.472136 | 4.472136 | 2, 6.708204 | 7.338205 | 8.303899 | 1.86 |
| 5demicube | 26 | 1.224745 | 1.380281 | 1.224745 | 10, 1.414214 | 1.414214 | 2.874418 | 2.08 |
| 6demicube | 44 | 1.000000 | 1.151037 | 1.000000 | 12, 1.414214 | 1.414214 | 2.418735 | 2.10 |
| 7demicube | 78 | 0.845154 | 1.021907 | 1.000000 | 14, 1.414214 | 1.414214 | 2.183926 | 2.14 |
| 8demicube | 144 | 0.738549 | 0.936467 | 1.000000 | 16, 1.414214 | 1.414214 | 2.189167 | 2.34 |
| 9demicube | 274 | 0.661438 | 0.873506 | 1.000000 | 18, 1.414214 | 1.414214 | 2.075661 | 2.38 |
| 10demicube | 532 | 0.603023 | 0.823943 | 1.000000 | 20, 1.414214 | 1.414214 | 2.044584 | 2.48 |

### 4C.8 The server, the work, the kernel and the search

- **Rule** (`ndServerRule`): P at most 24, the form (cell 0 to C - 1, u of d - 1 numbers at most 1e4 in absolute value, no `q`), s at least the problem's bound less 1e-9 relative, then `validateND` at tolerance 1e-9 with 20 reasons and the step budget; an entry on a problem without a record must also pass at 1e-11 (section 12).
- **Work budget.** A check counts steps: a node created, a projection solved, a join candidate examined. It stops after **ND_MAX_WORK = 60,000** steps and refuses the configuration with a plain reason ("the check reached its work limit ... too many balls crowd the (d-3)-faces"): sound, incomplete. [CALC, Apple Silicon, one core] A step costs 5 to 26 microseconds, so a refused check costs at most about a second and a half. Measured on the test packings (touching, tolerance 1e-9):

  | container | P = 1: steps / ms | layout: P: steps / ms | P = 4 | P = 24 |
  |---|---|---|---|---|
  | 5simplex | 1: 7 / 0.09 | 6: 12 / 0.10 | 4: 28 / 0.15 | 24: 78 / 0.21 |
  | 6simplex | 1: 8 / 0.02 | 7: 14 / 0.07 | 4: 32 / 0.07 | 24: 78 / 0.12 |
  | 7simplex | 1: 9 / 0.03 | 8: 16 / 0.03 | 4: 8 / 0.02 | 24: 76 / 0.12 |
  | 8simplex | 1: 10 / 0.04 | 9: 18 / 0.06 | 4: 8 / 0.02 | 24: 84 / 0.21 |
  | 9simplex | 1: 11 / 0.06 | 10: 20 / 0.08 | 4: 8 / 0.03 | 24: 78 / 0.23 |
  | 10simplex | 1: 12 / 0.08 | 11: 22 / 0.10 | 4: 8 / 0.03 | 24: 81 / 0.25 |
  | 5orthoplex | 1: 7 / 0.03 | 2: 4 / 0.02 | 4: 20 / 0.07 | 24: 81 / 0.27 |
  | 6orthoplex | 1: 8 / 0.01 | 2: 4 / 0.01 | 4: 17 / 0.02 | 24: 94 / 0.13 |
  | 7orthoplex | 1: 9 / 0.02 | 2: 4 / 0.01 | 4: 19 / 0.03 | 24: 90 / 0.15 |
  | 8orthoplex | 1: 10 / 0.03 | 2: 4 / 0.01 | 4: 16 / 0.03 | 24: 91 / 0.29 |
  | 9orthoplex | 1: 11 / 0.05 | 2: 4 / 0.01 | 4: 19 / 0.05 | 24: 90 / 0.22 |
  | 10orthoplex | 1: 12 / 0.06 | 2: 4 / 0.01 | 4: 20 / 0.05 | 24: 92 / 0.28 |
  | 5demicube | 1: 18 / 0.04 | 10: 20 / 0.06 | 4: 72 / 0.19 | 24: 76 / 0.18 |
  | 6demicube | 1: 988 / 2.8 | 12: 24 / 0.07 | 4: 2608 / 8.3 | 24: 115 / 0.39 |
  | 7demicube | 1: 1934 / 8.1 | 14: 28 / 0.13 | 4: 56 / 0.87 | 24: 226 / 1.1 |
  | 8demicube | 1: 5392 / 39.8 | 16: 32 / 0.24 | 4: 64 / 1.00 | 24: 119 / 0.92 |
  | 9demicube | 1: 18 / 0.59 | 18: 36 / 0.55 | 4: 72 / 2.4 | 24: 90 / 1.3 |
  | 10demicube | 1: 20 / 1.4 | 20: 40 / 1.2 | 4: 80 / 5.3 | 24: 97 / 2.7 |
  The slowest checks are balls at the centre of a demicube facet with the edge at 1 + 1e-9 (the P = 1 packings and the P = 4 on the 6-demicube): dozens of hinges at the same distance 1/2 s from the centre, so many chains survive the pruning (the 6-demicube: 988 and 2,608 steps; the 8-demicube: 5,392 steps, 41 ms). Seeded runs of the search gave at most 7,824 steps (6-demicube, P = 12) and 19 ms.
  **Ties are the cost.** A ball whose hinges are *exactly* 1/2 away is the worst case, and the centre of a demicube facet at s = 1 has dozens of them at that distance: the engine's reach is 1/2 (1 + 1e-12), so a hinge within 5e-13 of 1/2 ties, and every chain around it survives the pruning. One ball at the centre of a facet, steps (time) at s - 1 = 1e-10 / 1e-12 / 1e-13: the 6-demicube 988 (2 ms) / 68,948 (167 ms) / 68,948 (294 ms);
  the 7-demicube 1,934 (11 ms) / over 5,000,000 (25 s) / 204,110 (1.3 s); the 8-demicube 5,392 (38 ms) / over 5,000,000 (66 s) / over 5,000,000 (58 s). With the cap of 60,000 steps the server **refuses** such a ball for its work limit (sound, incomplete); a margin of 1e-10 and more is decided at once. The reach was therefore kept below the 1e-12 margins of the packings the server is to accept.
- **Heavy queue.** `boundaryIsHeavy` estimates a check at ND_BALL_MS per ball: three times the largest cost per ball measured (0.08 ms on the simplex and orthoplex, 5.5 on the 6-demicube, 9.5 on the 7-, 41 on the 8-, 0.6 on the 9-, 1.4 on the 10-demicube), rounded up: 0.3, 17, 30, 125, 2 and 4 ms (the 5-demicube 0.3). A check of 50 ms or more waits in the server's heavy queue, one at a time, and is charged its estimate:
  **the 6-demicube from P = 3, the 7-demicube from P = 2, the 8-demicube always, the 10-demicube from P = 13**; nothing else.
- **Kernel** (section 14). `boundaryKernel` returns `ndKernel` for these containers: n variables per piece, the energy of `nd/model.ts` (squared violations of the clearance, of every pair and loop below 1, found by the engine, with an analytic gradient in the displacements through the developments' linear parts; the size derivative is a central difference at fixed local coordinates plus the centres' own motion), weights, focus and the search's
  finish (`finalizeState` with shrink). A ball nearer than 1/2 to a hinge of its own cell is not developed (its pair terms use the chord in R^D, a lower bound). Seeds: the exact layout where one is derived (P = 1; simplex P <= d + 1; demicube P <= 2d; orthoplex P = 2) and a sequential random layout otherwise (balls added at random points, accepted when they conflict with none, the edge raised 6% after a failed round). Random starts put the balls on distinct cells at 1.1 times the layout edge, the demicube's balls on its facets, not its simplex cells.
  The kernel contract (`checkKernel`) passes on 8 containers spread over the three families and d = 5 ... 10; an energy evaluation with its gradient costs 0.02 to 0.1 ms per ball on the simplex and orthoplex, 0.05 to 1.7 ms on the demicube (`BOUNDARY_ENERGY_MS` in `apps/web/src/space/model.ts`).
- **Search.** `searchBoundary` on these problems seeds and finishes (no optimiser: the problem sets ship without records); the generic methods run on the kernel for anything more.
- **Tighter** must not see these containers: the game's slices and id patterns do not name `boundary` (`apps/web/test/tighter-boundary-excluded.test.ts` pins that for the 432 ids as it does for the others).

### 4C.9 On the site

The Problem set menu's option of each of these dimensions is now "4D on 5D: Polytope boundaries" ... "9D on 10D: Polytope boundaries" (value `boundary-<d>cube`, unchanged, so old links keep working; their address and the Records page's tab still open the d-cube first, the family's `open` filter). A family's picker offers its four containers (the d-cube
first, then the simplex, orthoplex and demicube); balls only on the last three. The board draws a packing as a table of the balls (index, facet, centre in the facet, centre in x1 ... xD) and a chart of the centres projected onto two axes the user picks, each ball in the colour of its facet, inside the outline of the projected container (`apps/web/src/space/render/ndpoly.ts`); a ball in a demicube facet
x_i = const is a disk when neither axis is i.

### 4C.10 Limitations

- **Balls only.** Polytope pieces on these boundaries are a later phase: the engine's chains, windows and joins are for points; a polytope piece needs its parts clipped by the windows in d - 1 dimensions (the 3-dimensional `parts.ts` and `polyhedron.ts` do not go there) and an exact separating-axis test in Q(sqrt k).
- **No records**, test packings only (4C.7); the 24-ball packings are loose (1.7 to 2.5 times the bound), the optima are unknown.
- **The demicube's one-ball bound** is the proved one (4C.5), not the sharp one for d >= 7; the other bounds are the volume bound and the one-ball bound only.
- **Finiteness of the chain enumeration** is not proved (4C.2); the budget makes it safe. A configuration with many balls crowded against hinges is refused with a reason rather than decided (sound, incomplete); packings with margins of 1e-12 and more are decided.
- **The exact certifier is slow on the big demicubes** (seconds to tens of seconds per check at d = 7 ... 10), because its numbers are Python integers in Q(sqrt k); the engine takes milliseconds.
- **The vertex rule** (section 3) is argued, not proved, in every d.
- Not covered: the Gosset polytopes (2_21, 1_22, 3_21, 4_21) of dimensions 6 to 8, and the 4-polytopes' rotated Platonic pieces (4A.8).


## 5. Lower bounds

Only proven bounds are stored [CALC]. Each is stored as a binary64 value rounded down, and the
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
- The d-simplex, d-orthoplex and d-demicube (d = 5 to 10, section 4C) take balls only; their problems have lower bounds and test packings, no records, and the demicube's one-ball bound is not known to be sharp for d >= 7.
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
- The d-cubes (section 4B): `packages/boundary/src/ncube/` (`lattice.ts` words and the reach search, `validate.ts`, `bounds.ts`, `move.ts`, `model.ts`, `kernel.ts`, `seeds.ts`, `rule.ts`, `search.ts`),
  `python/boundary/ncube.py` (the exact certifier), `python/boundary/test_ncube.py` and `packages/boundary/test/ncube.test.ts` (tests), `python/boundary/build_ncube_packings.py` (the test packings),
  `packages/boundary/scripts/gen-ncube-crosscheck.ts` and `gen-ncube-bounds.ts` (the fixtures `ncube_crosscheck.json`, `ncube_bounds.json`).
- The simplex, orthoplex and demicube (section 4C): `packages/boundary/src/nd/` (`names.ts`, `container.ts` the cell complexes, `chains.ts` half-chains, windows, joins and the hinge distance, `move.ts` the walks, `validate.ts`, `bounds.ts`, `model.ts`, `kernel.ts`, `seeds.ts`, `rule.ts`, `search.ts`),
  `python/boundary/ndpoly.py` (the exact certifier), `python/boundary/test_ndpoly.py` and `packages/boundary/test/nd.test.ts` (tests), `python/boundary/build_nd_packings.py` (certifies the test packings),
  `packages/boundary/scripts/gen-nd-crosscheck.ts`, `gen-nd-bounds.ts` and `gen-nd-packings.ts` (the fixtures `nd_crosscheck.json`, `nd_bounds.json`, `nd_packings_raw.json`), `apps/web/src/space/render/ndpoly.ts` (the board's chart and table), `apps/server/test/spaces-boundary-nd.test.ts`.

```
python3 python/boundary/certify.py data/boundary/records.json
python3 python/boundary/test_boundary.py
python3 python/boundary/test_general.py
python3 python/boundary/test_ncube.py
python3 python/boundary/test_ndpoly.py
npx tsx packages/boundary/scripts/search-sources.ts
python3 python/boundary/build_records.py
npx tsx packages/boundary/scripts/check-records.ts
npx tsx packages/boundary/scripts/gen-crosscheck.ts
npx tsx packages/boundary/scripts/gen-containers.ts
npx tsx packages/boundary/scripts/gen-general-crosscheck.ts
npx tsx packages/boundary/scripts/gen-general-packings.ts
python3 python/boundary/build_general_packings.py
python3 python/boundary/build_ncube_packings.py
npx tsx packages/boundary/scripts/gen-ncube-crosscheck.ts
npx tsx packages/boundary/scripts/gen-ncube-bounds.ts
npx tsx packages/boundary/scripts/gen-nd-crosscheck.ts
npx tsx packages/boundary/scripts/gen-nd-bounds.ts
npx tsx packages/boundary/scripts/gen-nd-packings.ts
python3 python/boundary/build_nd_packings.py
```

Timing (Apple Silicon, shared machine): the 11 searches about 2.5 minutes in all (the longest,
balls P = 13 and 15, 44 and 49 s); building the records file 3 s; certifying all 32 records with
their exact constructions 2.3 s; the Python tests about 26 s; the engine check 0.3 s for 64
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
- H. S. M. Coxeter, *Regular Polytopes*, 3rd ed. (Dover, 1973), 7.2 and 8.7: the simplex, the orthoplex (cross polytope) and the demicube (half-cube) in every dimension, their
  vertices and facets, used by section 4C (the facet counts are the book's; the combinatorics were also read from the vertex sets by brute force).

## 12. On the site: the API and its rule

The server lists the space as `boundary` (title "Polytope boundaries", viewer
`/spaces/boundary.html`) through the same routes as the other spaces; the API page, `/api`,
describes them all. Everything below goes through the spaces registry, `packages/spaces`.

- `GET /api/spaces`: the space with its conventions and its 152 problems (its `domain`: polytope, piece, P; the pieces and the P range depend on the polytope).
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

1. The id must be a problem of the domain, in its one written form (no leading zeros): balls and cubes on the tesseract, P = 1 to 16, balls on
   the 5-cell, 16-cell, 24-cell, 120-cell and 600-cell, P = 1 to 24, balls and cubes on the d-cubes (section 4B), balls on the d-simplex, d-orthoplex and d-demicube, d = 5 to 10, P = 1 to 24 (section 4C). P is at most 64 for the space and at most the
   catalogue's largest P (24) for a submission.
2. The configuration must be well formed, or the answer is 400 with what is wrong: s a finite
   number in (0, 1e4]; `pieces` a list; each piece `{cell, u}` with `cell` an integer from 0 to C - 1 (C = 8 on the tesseract, 5, 16, 24, 120, 600 on the others)
   and `u` three finite numbers, each at most 1e4 in absolute value; a cube may have `q`, four finite
   numbers whose norm is within 1e-9 of 1; a ball has no `q`. Unknown fields are dropped.
3. A wrong number of pieces, or an s below the problem's lower bound (less a relative 1e-9), is an
   invalid packing: 422 with the reason, before any geometry.
4. `validateBoundary` (section 4.5's engine validator) at the touching tolerance 1e-9, with at most
   1024 chains of cells per piece: 422 with up to 20 reasons when it refuses.
5. An entry below score_cert by more than a relative 1e-9 must also pass at tolerance 1e-11
   ("beats the research record only within the float tolerance" otherwise). Every record passes it
   at its stored s. On a problem without a record (the five new containers) every entry must pass it
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
that problem. The Problem section offers Polytope (Tesseract), Piece (Balls, Cubes) and Pieces (P),
each with the values the catalogue has. The stage shows the drawing (below), the panel the usual
values (the edge s of the packing shown, the certified value, the status, the lower bound, the
site best, whether the problem takes entries) and two of this kind:

- **Without folding**: the nofold value of section 6, with its derivation under "Why".
- **Folding gains**: 1 - best / nofold, where best is the problem's best edge (the research record
  or a better site entry), to one decimal of a percent. For a record this is section 7's gain, to
  within the 1e-12 of certification; a site entry that beats the record raises it.

Run optimizer runs section 12's search in the browser, one complete search call per round, of
300,000 energy evaluations (up to three times that while rounds get nowhere), about 1.5 to 2.2 s
each at P = 16, and Submit sends a result that passes the server's rule
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
