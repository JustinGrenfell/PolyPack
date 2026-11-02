# 3D on 4D, polytope boundaries: balls (diameter 1) and cubes on the boundary of a regular 4-polytope

**Dimension.** Under Dimensions this is **3D on 4D**: solid 3D pieces packed on the 3-dimensional
boundary of a 4D polytope, with no container around them. Problem ids and records use the space
name `boundary`.

P balls of diameter 1 (radius 1/2) or P unit cubes (edge 1) are packed on the boundary of a tesseract of
edge s. Since 2026-11-01 balls can also be packed on the boundary of the other five regular
4-polytopes, the 5-cell, 16-cell, 24-cell, 120-cell and 600-cell (section 4A). The score is the smallest s,
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
within about 1e-16 s of the engine's. A centre that falls outside its closed cell by less than 1e-12 of an edge is moved onto the facet planes it
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
(lam = 1 + k 1e-12; k = 0 or 1 for all) and certified strictly valid by `general.py`. They are not good packings (no search ran) and no optimality is claimed.

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
- The ball centre is read through float frames (4A.3): the certified configuration is within about 1e-16 s of the one the engine sees.
- The 120-cell and 600-cell are nearly flat. Their problems are open, and their bounds are those of 4A.4; nothing is known about their optimal packings.

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

```
python3 python/boundary/certify.py data/boundary/records.json
python3 python/boundary/test_boundary.py
python3 python/boundary/test_general.py
npx tsx packages/boundary/scripts/search-sources.ts
python3 python/boundary/build_records.py
npx tsx packages/boundary/scripts/check-records.ts
npx tsx packages/boundary/scripts/gen-crosscheck.ts
npx tsx packages/boundary/scripts/gen-containers.ts
npx tsx packages/boundary/scripts/gen-general-crosscheck.ts
npx tsx packages/boundary/scripts/gen-general-packings.ts
python3 python/boundary/build_general_packings.py
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
   the 5-cell, 16-cell, 24-cell, 120-cell and 600-cell, P = 1 to 24. P is at most 64 for the space and at most the
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
