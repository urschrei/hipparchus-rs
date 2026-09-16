# Design

## 1. Purpose

This crate reimplements the Hipparchus geopositioning model (Lukatela 1987;
Lukatela & Russell 1992) from its published descriptions, so that the ideas
can be studied, tested and extended. This document records the model as we
understand it, the design of v0.1, what we verified and how, the problems we
found in the 1992 code, and a roadmap detailed enough to resume work without
re-reading the sources.

The [README](README.md#references) lists the references. "AC8" below means
the Auto-Carto 8 paper, "DDJ" the 1992 *Dr. Dobb's Journal* article and its
companion archive, and "DGG" the 2000 terrain-model paper.

## 2. Background

### 2.1 The model

AC8 describes three components:

1. **Global coordinates as direction cosines** of the ellipsoid normal. The
   domain is restricted to the surface, there are no singular points, and
   64-bit floats give sub-millimetre resolution. Geodetic problems are solved
   by iteration rather than by series expansion in the eccentricity.
2. **A spheroidal Voronoi cell structure.** A cell is the set of surface
   points "closer" to its centre than to any other centre (the quotation
   marks are in the original, and matter; see R5). The structure is stored
   as the centres' global coordinates plus a circular list of neighbour
   identifiers per cell, i.e. the Delaunay dual. Centres are placed by the
   application: by data density, activity, size limits, or existing control
   points.
3. **Cell-local coordinates.** A point is a cell identifier plus offsets in
   an "ellipsoid-modified stereographic" plane tangent at the cell centre.
   AC8 gives 80 bits per point (16-bit cell, 32-bit local values) at
   sub-millimetre resolution. Local coordinates are for compact storage and
   display only; metric problems use global coordinates.

Objects (AC8, DDJ):

- **Points:** cell identifier plus local coordinates.
- **Lines:** ordered list of traversed cells, each with its vertices;
  boundary-crossing points may be stored permanently.
- **Regions:** directed circular boundary plus an encoded list of cells
  wholly inside. Many relationships are then decided from cell lists alone.
  On a sphere the complement of a region is finite.

Operations:

- **Home-cell search:** start from an application-chosen cell; move to the
  nearest neighbour centre while one is nearer than the current centre.
- **Restructuring:** a cell that becomes too dense can be split, with changes
  confined to its neighbourhood.
- DGG adds a seamless global TIN built cell by cell, with interpolated
  elevations on cell edges.

### 2.2 Available resources

- The library (`hipparch.lib`) was distributed as a binary. The auxiliary
  library and utilities, including the cell-structure builder, were C source
  but only within the commercial kit.
- The Geodyssey archive site states that exclusive rights to the libraries
  were sold to Microsoft and that the company was voluntarily dissolved.
- Public code: the DDJ companion archive only (`ALGEBRA.C`, `BUILD.C`,
  `VIEW.C`, coastline data). It is the spherical, object-centred precursor;
  it does not contain the Voronoi engine.
- Public documentation: AC8, DDJ, DGG, the tutorial, the function reference
  and the utilities guide (links in the README). The function reference is
  the best source for the real library's data structures.

## 3. Crate overview (v0.1)

The README lists the modules.

Dependencies: `geo` (`Point` conversions and Delaunay via spade; Haversine
in tests), `spade` (the `Triangulation` trait and vertex handles for the
triangulation that `geo` returns; the requirement must stay compatible with
`geo`'s), `robust` (exact `orient3d`).

Conventions:

- Coordinates are geocentric on a unit sphere; lat/lon are degrees.
  `EARTH_RADIUS` is 6,371,008.8 m, matching `geo`'s Haversine.
- Triangles are counter-clockwise seen from outside the sphere.
- Cell identifiers are `u32` indices into the centre list.

## 4. Core representation

### 4.1 Direction cosines

`V3::from_lat_lon` treats latitude as geocentric, as `ALGEBRA.C` does.
Nothing in v0.1 depends on the ellipsoid. `from_lat_lon` and `normalized`
return unit vectors; `V3::new` and the arithmetic operators do not
normalise.

### 4.2 Distance

- `arc(a, b) = atan2(|a × b|, a · b)`. Well conditioned at all separations,
  unlike `acos(a · b)`. Verified against `geo`'s Haversine to 1 µm.
- `arc_ddj` reproduces `ArcDist`: `c + c³/24` with `c` the chord. These are
  the first two terms of `2 asin(c/2)`. It underestimates long arcs (2.33
  against π at the antipode). DDJ used it only for culling, where an
  underestimate is safe.
- The home-cell walk and cell assignment use the dot product, which is
  monotone in arc length and needs no square root or trigonometry.

### 4.3 Stereographic frame

For a centre `c`:

- east `e = (−c_y, c_x, 0) / hypot(c_x, c_y)`; north `n = c × e`.
- `e × n = c`, so (`e`, `n`, `c`) is right-handed with `c` outward:
  counter-clockwise in the plane is counter-clockwise seen from outside.
- `map(p) = 2 / (1 + c·p) · (e·p, n·p)`: projection from `−c` onto the plane
  tangent at `c`. Singular only at `p = −c`.
- `unmap(E, N) = ((4 − ρ²) c + 4 (E e + N n)) / (4 + ρ²)`, with
  `ρ² = E² + N²`. The implementation normalises the result, as the original
  does.
- A point at arc `r` from `c` maps to plane distance `2 tan(r/2)`
  (`plane_extent`), which is at least `r`. The linear scale factor is
  `sec²(r/2) ≥ 1`, so one plane unit is at most one unit of ground distance.

This is algebraically the same mapping as `MapStereo`/`UnMapStereo`; the
original constructs it geometrically (line–plane and line–sphere
intersections). The tests confirm agreement to 1e-12.

**Pole handling.** The original divides by `hypot(c_x, c_y)`, which is zero
at the poles. When that value is at most 1e-300 we use `e = (0, 1, 0)`. The
frame is ill-conditioned near a pole (a small change in `c` rotates it), but
it is a deterministic function of the stored centre, so encoding and
decoding against the same centre always agree. `poles_round_trip` covers
both poles.

## 5. Cell structure

### 5.1 Definition and invariants

`CellStructure::new(centres)` holds, for each cell `i`:

- `centres[i]`: the site.
- `neighbours[i]`: Delaunay neighbours in counter-clockwise order.
- `vertices[i]`: Voronoi vertices in the same order; vertex `k` lies between
  neighbours `k` and `k + 1`.
- `radius[i]`: maximum arc from the centre to any point of the cell, which
  for a spherically convex cell is attained at a vertex.
- `frames[i]`, `scale[i]`: encoding parameters (5.7).

Invariants:

- I1. Every centre is a Delaunay vertex (no duplicates, no `TooClose`).
- I2. Each neighbour ring is closed and consistent: `j ∈ N(i) ⇔ i ∈ N(j)`.
- I3. Every point `p` with home cell `i` satisfies `arc(p, c_i) ≤ radius[i]`
  up to rounding (tested with a 1e-12 relative margin).
- I4. `radius[i] ≤ MAX_CELL_RADIUS` (2.0 rad).
- I5. The structure is immutable; any change means a rebuild (see R4).

### 5.2 Delaunay triangulation

The spherical Delaunay triangulation of unit vectors is their convex hull.
Two backends produce it.

**`sphere_delaunay_geo` (default).**

1. Project every site except site 0 with `Frame::new(−s₀)`: tangent at the
   antipode, projection from `s₀`.
2. Stereographic projection maps circles not through `s₀` to circles and
   circles through `s₀` to lines. An empty spherical cap bounded by the
   circle through `a, b, c` therefore corresponds to an empty planar disk,
   so the planar Delaunay triangles are the spherical ones not incident to
   `s₀`.
3. Triangles incident to `s₀` correspond to lines through two projected
   sites with all others on one side: the planar convex hull edges. spade
   returns hull edges clockwise; each edge `(from, to)` gives the triangle
   `(from, to, s₀)`.
4. The frame at `−s₀` is right-handed seen from outside and the projection
   is conformal and orientation-preserving away from `s₀`, so spade's
   counter-clockwise faces need no correction.
5. Sites are inserted in Morton order (21 bits per axis of the unit vector)
   because `geo`'s `unconstrained_triangulation_raw` inserts sequentially
   and spade starts each point location from the previous insertion. This
   took 10⁶ random sites from 77 s to 0.69 s.
6. The raw variant is used deliberately. The configured variant snaps
   coordinates within 1e-4, which on a unit-sphere plane merges sites about
   600 m apart.

Errors:

- `Duplicate(i)`: site `i` is the first site equal to an earlier site. Both
  backends report the same index.
- `TooClose(i)`: two distinct sites project to identical coordinates; `i`
  is the lowest index that coincides with an earlier site.
- `Projection` when spade rejects a coordinate (magnitude above 2²⁰¹, or
  not finite). Sites near `s₀` project to about `4/θ`, so the limit is not a
  practical concern. spade also rejects non-zero magnitudes below 2⁻¹⁴²;
  such coordinates become zero, a change below the rounding error of
  `map`.
- `TooFewSites` for fewer than four sites or all projected sites collinear.

Exactness: spade's predicates are exact on the projected floating-point
coordinates, not on the original vectors, so near-cospherical ties can
resolve differently from an exact 3D test. The sites +X, +Z, −Z, +Y and
`(−1/3, 2/3, 2/3)` show this: the last site lies within rounding of the
circle through +X, +Z and +Y, and the two backends choose different
diagonals. The result is valid to within rounding (every `orient3d` is
above −1e-14). See O6.

**`sphere_delaunay` (reference).** Incremental 3D hull:

- Initial tetrahedron from sites 0, `i1`, `i2` (three distinct sphere points
  are never collinear) and the first `i3` off their plane; orientation fixed
  with one `orient3d`.
- Face `(a, b, c)` is visible from `p` if `orient3d(a, b, c, p) < 0` (robust's
  sign convention: positive means below a counter-clockwise face).
- Horizon: visible directed edges whose twin is not visible. Each gives the
  new face `(a, b, p)`, which keeps the old orientation.
- Coplanar faces (`orient3d == 0`) count as not visible. Visibility is
  constant over coplanar faces, so the visible set is always a disc and the
  result remains a valid triangulation of the tie.
- A site with no visible face is a `Duplicate` or, if distinct, `TooClose`:
  it lies within rounding of the hull, which requires a separation below
  roughly a metre on the Earth.
- Cost is O(n·F), i.e. O(n²): 6.7 s at 20,000 sites. It remains useful as an
  exact oracle and as the basis for local restructuring (R4).

### 5.3 Neighbour rings and Voronoi vertices

For each triangle `(a, b, c)` the circumcentre is
`normalize((b − a) × (c − a))`. With counter-clockwise faces this normal
points outward, so it is the Voronoi vertex on the correct side even for
very large triangles. Site `a` collects `(b, c, circumcentre)` from each
incident triangle `(a, b, c)`; the ring follows `c` to the triangle whose
`b` equals it. Each step searches the fan, so a ring costs the square of
the degree, which is small.

### 5.4 Cell radius

`radius[i]` is the largest arc from the centre to its Voronoi vertices.
`MAX_CELL_RADIUS = 2.0` rejects cells so large that the stereographic scale
becomes unusable (`plane_extent` grows without bound towards π). Global
seeds prevent this in practice. If gnomonic per-cell geometry is adopted
(R3), the limit must drop below π/2.

### 5.5 Home-cell walk

```
cur = start; best = p · c[cur]
loop: move to the neighbour with the largest p · c[nb] if it exceeds best
```

Correctness: the Voronoi cell of `v` is the intersection over its Delaunay
neighbours `u` of the hemispheres `{x : x·v ≥ x·u}`. If `p` is not in `v`'s
cell, some Delaunay neighbour is strictly nearer, so the walk only stops at
the nearest centre. Termination follows from the strictly increasing
`best`. `walk_finds_nearest_centre` checks this against brute force on
5,000 queries.

Ties: a point equidistant from two centres stays wherever the walk reaches
first, so its home cell can depend on `start`. This is harmless for
encoding (both cells decode correctly) but would matter for a unique stored
representation and for deduplication. See open issue O1.

Hints: the caller supplies `start`. For sequential vertices the previous
vertex's cell is almost always correct (0.03 steps per point on the
coastline data, against 7.7 from a fixed start).

### 5.6 Cap queries

`cells_in_cap(q, r, hint)` walks to `q`'s home cell, then searches outward,
keeping a cell if `arc(q, c) ≤ r + radius(c)` and expanding only from kept
cells.

- The test is necessary for a cell to meet the cap (triangle inequality),
  so every intersecting cell passes it.
- The cells meeting a cap form an edge-connected set: the cap is connected
  and the cells tile the sphere. Since that set is contained in the passing
  set, the search reaches all of it.
- False positives are possible and harmless: a cell can pass without
  meeting the cap.

The Dublin 500 km query returned 9 cells and read 2.9% of the coastline
vertices, finding all 662 true hits.

### 5.7 Local encoding

- `scale = L / plane_extent(radius)` with `L = 2 × 10⁹`, leaving headroom
  below `i32::MAX`.
- A point in its home cell has arc at most `radius` (I3), so each plane
  coordinate has magnitude at most `plane_extent(radius)` and each scaled
  value at most `L`. Overflow is impossible by construction.
- Values are rounded to nearest, not truncated.
- Plane resolution is `plane_extent(R) / L` radians per unit, and ground
  resolution is no worse (4.3). Rounding contributes at most half a unit per
  axis, so the worst position error is `plane_extent(R) · √2 / (2L)`. For
  the largest coastline cell (R = 1,690 km) the bound is 0.60 mm; the
  measured worst case was 0.58 mm.
- `Local` is 12 bytes (`u32` cell, two `i32`). AC8's layout used a 16-bit
  cell identifier (10 bytes).

## 6. Adaptive builder

`build::adaptive(points, capacity, seeds)`:

1. Start from `seeds` Fibonacci-lattice centres, so empty regions are still
   covered by cells of bounded size.
2. Build the structure and assign every point with a walk hinted by the
   previous point's cell.
3. Replace each cell with more than `capacity` points by the two centres of
   a spherical 2-means split (initialised from the point farthest from the
   mean and the point farthest from that; 8 iterations).
4. Repeat until no cell needs splitting or no overloaded cell can be split
   (its points are identical, or too close for dot products to separate;
   see O8).

Limitations:

- Every round rebuilds the whole structure; AC8 emphasises local change.
- Splitting replaces a centre, which requires deletion in any incremental
  scheme (R4).
- Cells are never merged, and nothing prevents two new centres from landing
  within a metre of an existing one (which would return `TooClose`).
- Point sets are held in memory per round.

## 7. The 1992 code

### 7.1 File formats

- `COAST*.IDX`: 40-byte little-endian records from the 16-bit MS C 6.0
  compiler: centre (3 × `f64`), radius in radians (`f64`), offset into
  `COAST.DAT` (`i32`), vertex count (`i16`), segment id (`i16`). A modern
  `struct` has a different layout, so both `ddj.rs` and `ref.c` read fields
  individually.
- `COAST.DAT`: `(i16, i16)` pairs, `value = (short)(plane × 32000 / radius)`.
  The sparse set is written first, then the dense set; both index files
  point into the same file. The indexes list 63,954 vertices in total,
  exactly matching `COAST.DAT`'s size.
- `COAST*.LNS`: ASCII; `*` lines start segments; one `lat lon` pair (degrees)
  per line. `COAST0.LNS` vertices correspond one to one, in order, with the
  dense index.

### 7.2 Findings

- **Precision.** The 16-bit encoding's worst error against the ASCII source
  is 111 m. The scale is per object (`32000 / radius`), and the cast
  truncates, so the error per axis is up to one unit: for the largest dense
  segment (radius 0.4096 rad, 2,609 km) one unit is 82 m and the bound is
  115 m. The median dense segment (39 km) has a unit of about 1.2 m.
- **Overflow.** The plane offset is `2 tan(r/2)`, not `r`, so values exceed
  32,767 once the object radius passes about 0.529 rad (3,370 km). The
  supplied data peaks at 32,455. Nothing checks this, and `ArcDist`'s
  underestimate of the radius makes it slightly worse.
- **Poles.** `MapStereo` and `UnMapStereo` divide by `hypot(x, y)`; a
  centre at a pole produces NaNs. The file header says "fuzz control" was
  removed from the library code, which presumably covered this.
- **Culling.** `VIEW.C` scans the whole index on every redraw. DDJ presents
  this as the problem the Voronoi structure solves.
- **Double projection.** Each vertex is unprojected from its object centre
  and reprojected to the display centre.
- **Sphere.** `LatLongToDcos3` is geocentric; the ellipsoidal algebra is only
  described in prose.

### 7.3 Reference data

`tools/ddj-reference/ref.c` calls the original `ALGEBRA.C` (not included;
fetch instructions in that directory's README) and writes:

- `tests/data/stereo.csv`: 2,000 cases of centre, point, `MapStereo`
  result, `UnMapStereo` of that result, and `ArcDist`. Centres are within
  89° of the equator (the original fails at the poles); points are up to
  166° from the centre (mean 42°).
- `tests/data/coast0_first.csv`: the first decoded vertex of every dense
  segment, decoded as `VIEW.C` does.

## 8. Verification

Tests marked (P) are property tests using Hegel. Their unit vectors come
from arbitrary coordinates in [−1, 1], so the inputs include poles, axis
points and cocircular sets. Tests that build a cell structure use a
Fibonacci lattice plus drawn sites, and drop any site within 2e-7 rad of
another (the `TooClose` precondition). Set `HEGEL_TEST_CASES` to run more
cases than the default 100.

| Test | Establishes |
|---|---|
| `stereo_matches_original_c` | `Frame::map`, `Frame::unmap` and `arc_ddj` agree with the C to 1e-12 (1e-14 for `arc_ddj`). |
| `coast_decode_matches_original_c` | The index reader reproduces the C offsets, and the decoder reproduces the first vertex of each of the 992 dense segments to 1e-12. |
| `poles_round_trip` | The pole frame round-trips at both poles. |
| `lat_lon_round_trips` (P) | `lat_lon` and `from_lat_lon` round-trip to 1e-14. |
| `arc_agrees_with_geo_haversine` (P) | `arc` matches `geo` to 1 µm. Near the antipode the tolerance grows with Haversine's conditioning. |
| `arc_ddj_does_not_exceed_arc` (P) | `arc_ddj` does not overestimate, so it is safe for culling. |
| `frame_round_trips` (P) | `unmap(map(p))` returns `p` to 1e-14 for arcs up to `MAX_CELL_RADIUS`. |
| `frame_maps_arc_to_plane_extent` (P) | `map` puts a point at plane distance `plane_extent(arc)`, as the overflow argument of 5.7 requires. |
| `clustered_sites_triangulate` | Both backends triangulate a 1 km cluster that contains the projection pole. |
| `delaunay_backends_are_valid` (P) | Each result is a closed surface with 2n − 4 faces and passes the empty-circumcircle check: exactly for `sphere_delaunay`, to 1e-14 for `sphere_delaunay_geo`. |
| `delaunay_backends_agree_without_near_ties` (P) | The backends give the same edges when no edge is within rounding of a tie. |
| `backends_report_the_first_duplicate` (P) | Both backends report the same `Duplicate` index. |
| `cospherical_and_duplicate_sites` | Icosahedron (five coplanar neighbours per vertex): valid triangulations and five-neighbour rings; duplicates reported by index. |
| `neighbour_rings_are_symmetric` (P) | I2. |
| `voronoi_vertices_are_nearest_to_their_sites` (P) | Vertex `k` is equidistant from the centre and neighbours `k` and `k + 1`, and no centre is nearer. |
| `walk_finds_nearest_centre` (P) | The walk reaches a nearest centre from any start. |
| `cell_radius_bounds_its_points` (P) | I3. |
| `encoding_error_is_bounded` (P) | The round-trip error is at most `√2/2` resolution units (5.7). |
| `encoding_does_not_saturate` (P) | Local values stay below `i32::MAX` in magnitude. |
| `cap_query_contains_every_home_cell` (P) | A cap query contains the home cell of every point in the cap. |
| `adaptive_cells_respect_capacity` (P) | Loads sum to the point count, and only a cell with one point exceeds the capacity. |
| README doctest | The README usage example compiles and round-trips below 1 mm. |

CI (`.github/workflows/ci.yml`) checks formatting, lints, tests and the
minimum supported Rust version (1.88, set by `geo`).

## 9. Performance

Release build on a shared Linux sandbox; indicative only.

| Task | Result |
|---|---|
| Delaunay, 1,000 random sites | exact 17 ms; `geo` 0.6 ms |
| Delaunay, 20,000 random sites | exact 6.7 s; `geo` 12 ms |
| Delaunay, 10⁶ random sites | `geo` 0.69 s (77 s without Morton ordering) |
| Adaptive build, coastline (52,766 points, capacity 256) | 26 ms |
| Encode and decode 52,766 points | 7 ms |

`examples/hull_bench.rs` and `examples/coast.rs` reproduce these.

## 10. Open issues

- O1. Tie-breaking in the walk is start-dependent. A unique stored
  representation needs a rule, e.g. among equidistant centres choose the
  lowest identifier, which requires an exact dot-product comparison (O3).
- O2. `encode` always walks from its hint, so a stale hint cannot cause
  overflow. Rounding can move a decoded point near a cell boundary into a
  neighbouring cell, so re-encoding it can change its cell.
- O3. There is no exact predicate for "`p` is nearer to `c` than to `u`",
  i.e. `sign(p · (c − u))`. `robust` provides only orientation and
  in-circle tests. Lines and regions (R1, R2) need it; Shewchuk-style
  expansion arithmetic would provide it.
- O4. `CellStructure` holds `Vec<Vec<_>>` rings. A compressed layout (offsets
  plus one flat array) would suit memory mapping (R6).
- O5. `HullError` and `BuildError` do not implement `std::error::Error` or
  `Display`.
- O6. `sphere_delaunay_geo` is Delaunay only to within rounding (5.2). A pass
  of exact `orient3d` edge flips after spade would make it exact.
- O7. `sphere_delaunay_geo` reports `TooClose` only when two sites project to
  identical coordinates. Distinct sites 1e-16 rad apart build a
  `CellStructure` whose Voronoi vertices are wrong. A minimum-length check on
  Delaunay edges would enforce the documented separation, because each
  site's nearest neighbour is a Delaunay neighbour.
- O8. Dot products do not separate points closer than about 1.5e-8 rad
  (10 cm), because `cos` is flat near zero. The walk and the adaptive split
  share this limit: `adaptive` keeps two points 1.5e-8 rad apart in one
  cell. The exact predicate of O3 would remove it.

## 11. Roadmap

### R1. Line objects (next)

Representation, following AC8:

```rust
struct LineSet { parts: Vec<Vec<CellRun>> }
struct CellRun {
    cell: u32,
    entry: Option<Crossing>,   // None at the start of a part
    vertices: Vec<(i32, i32)>, // local coordinates in `cell`
    exit: Option<Crossing>,    // None at the end of a part
}
struct Crossing { neighbour: u32, t: u32 } // fixed-point position on the shared Voronoi edge
```

The crossing is stored once per shared edge in a fixed direction (lower
cell identifier first) so that both cells agree exactly.

Traversal of a minor arc `a → b` starting in `a`'s home cell `c`:

1. The arc's great-circle normal is `g = a × b`. For each neighbour `u`, the
   bisector has normal `w = c − u`; points nearer `c` satisfy `x · w > 0`.
2. The arc leaves `c` through the bisector it crosses first. The crossing
   direction is `±normalize(g × w)`; keep the sign that lies on the arc,
   i.e. `(a × q)·g ≥ 0` and `(q × b)·g ≥ 0`.
3. Because the cell is an intersection of hemispheres and `a` satisfies all
   of them, the first bisector crossed along the arc is the exit edge.
4. Continue from the crossing in `u` until the current cell is `b`'s home
   cell.

Degeneracies to settle:

- The arc passes exactly through a Voronoi vertex (two bisectors tie):
  choose the lower identifier; the next step then leaves with zero length.
- The arc lies along a Voronoi edge: assign it to the lower identifier.
- A vertex lies exactly on an edge: covered by O1.
- Arcs of length π or more are ambiguous; require segments shorter than a
  fixed limit (for example π/2) and densify on input.

Predicates:

- Side of great circle `ab` for an input point: `orient3d(origin, a, b, x)`,
  exact via `robust`.
- Bisector side for an input point: exact dot sign (O3).
- Ordering of constructed crossings along the arc involves higher-degree
  polynomials. Plan: floating-point evaluation with a static error bound,
  falling back to exact expansion arithmetic.

Tests: traversals equal to a dense-sampling oracle (sample the arc finely
and take home cells); every stored crossing decodes to a point equidistant
from both centres; round trip of the coastline data as lines.

### R2. Regions

Representation: directed boundary rings with the interior on the left seen
from outside, plus a per-cell classification (`Inside`, `Outside`,
`Boundary`) stored as sorted identifier runs or a bitset. For each boundary
cell, store whether its centre is inside.

Construction:

1. Trace every ring with R1; the traversed cells are `Boundary`.
2. Determine one boundary cell centre's status from the orientation
   convention (the side of the first boundary segment).
3. Propagate status to other centres by crossing parity along arcs between
   centres. A Delaunay edge does not necessarily stay within the two cells it
   joins (obtuse triangles), so the parity count must use the boundary
   segments of every cell the arc traverses, which R1 provides.
4. Non-boundary cells take their centre's status.

Point in region: walk to the home cell. `Inside` or `Outside` answers
directly. For `Boundary`, count crossings of the arc from the cell centre
(known status) to the point with that cell's boundary segments. The arc stays
in the cell because cells are spherically convex and the two points are not
antipodal.

Validation: rings must be simple and must not cross one another;
multi-part regions and holes follow from the orientation convention. The
complement is formed by reversing every ring and swapping `Inside` and
`Outside`.

### R3. Set operations

Cell-level classification first:

| A | B | A ∩ B | A ∪ B |
|---|---|---|---|
| Outside | any | Outside | B |
| Inside | Inside | Inside | Inside |
| Inside | Boundary | B's geometry | Inside |
| Boundary | Boundary | geometry needed | geometry needed |

Geometry is only needed in cells that are `Boundary` in both, as in AC8
Figure 3. The measure worth reporting is the fraction of cells resolved
without coordinates.

For the per-cell geometry, stereographic coordinates are unsuitable because
great-circle arcs become circular arcs. A gnomonic projection about the cell
centre maps great circles to straight lines and cell edges to a convex
polygon, so `geo`'s boolean operations could clip within each cell. This
requires `MAX_CELL_RADIUS < π/2`. Rounding at cell edges must be controlled,
probably by using the stored `Crossing` points as shared exact vertices.

### R4. Local restructuring

- Insertion: walk to the new site's home cell, find a hull face visible from
  the new site among the triangles around that centre, grow the visible set
  by breadth-first search (it is a disc), and replace it with the horizon
  fan. This is the exact backend's step without the O(n) scan.
- Affected cells: the new cell and the horizon sites. Only points in the
  horizon cells need reassignment and re-encoding.
- Deletion (needed if splitting replaces a centre): re-triangulate the
  removed site's star polygon. Alternatively, split by keeping the old
  centre and inserting one new centre (insert-only).
- Identifiers: new cells get fresh identifiers; removed identifiers are
  retired, not reused, so stored `Local` values stay unambiguous.

### R5. Ellipsoid

- Store geodetic normals as the global coordinates, as AC8 describes;
  conversion to and from Earth-centred Cartesian positions goes through the
  prime-vertical radius of curvature.
- Distances and azimuths: implement the iterative methods if they can be
  reconstructed from AC8; otherwise use Karney's algorithms. `geo`'s
  `Geodesic` (via `geographiclib-rs`) is available for verification.
- The cell definition needs a decision. Voronoi cells by normal dot product
  are a valid tessellation (they are the spherical cells of the Gauss map)
  but are not cells by geodesic distance. AC8 puts "closer" and "distance"
  in quotation marks and says efficiency is very sensitive to the choice,
  which suggests a proxy. The two definitions differ by an amount related
  to the flattening (about 1/298); this should be quantified before
  choosing, and the choice documented.
- The "ellipsoid-modified stereographic" transformation is not described.
  A conformal sphere (for example Gauss–Schreiber) followed by the existing
  frame is one candidate.

### R6. Storage

- A binary format for cell structures and objects, suitable for memory
  mapping (Geodyssey used memory-mapped "Hipparchus binary objects", `.hbo`).
- Store a hash of the centre list with every object so that local
  coordinates are never decoded against the wrong structure.
- Options: 16-bit cell identifiers where fewer than 65,536 cells exist
  (10-byte points); delta and variable-length coding of vertices within a
  cell run.

### R7. Upstream to `geo`

`geo` has no spherical Delaunay or Voronoi. `sphere_delaunay_geo` is small,
uses only `geo` and spade, and has an exact reference test. A possible API:
a trait on `MultiPoint<f64>` (longitude, latitude) returning triangles as
index triples, with Voronoi cells as `Polygon`s. Cells crossing the
antimeridian or containing a pole need a policy before cells can be
returned as lon/lat polygons. The cell structure itself is better kept in
this crate.

### R8. Sources and contacts

- The VCF forum user `spd` (threads "A software thing I built: GPS on a
  25MHz 486-SX" and "Hipparchus GIS library docs / software from Geodesy ?")
  is looking for the same material.
- Hrvoje Lukatela hosts papers at lukatela.com and could clarify boundary
  handling, restructuring, the distance proxy and the ellipsoidal
  projection.
- Former customers' archives are the likeliest source of the auxiliary
  library and utilities source.

## 12. Decision log

| ID | Decision | Reason |
|---|---|---|
| D1 | `u32` cell identifiers | Removes the 65,536-cell ceiling; 16-bit remains an option for storage (R6). |
| D2 | Round local coordinates | Halves the worst quantisation error relative to the 1992 truncation. |
| D3 | `i32` locals with `L = 2 × 10⁹`, scale from `2 tan(R/2)` | Overflow-free by construction; worst error below 1 mm for cells up to about 2,800 km radius. |
| D4 | Fixed pole frame (east = +Y) | Deterministic; avoids the original's division by zero. |
| D5 | Dot product for nearest-centre decisions | Monotone in arc length; cheap; exact predicate still needed (O3). |
| D6 | `geo` Delaunay backend by default, exact hull as reference | Over 500× faster at 20,000 sites; identical results in tests except near-cospherical ties (O6). |
| D7 | Morton-ordered insertion | `geo`'s sequential insertion depends on input locality. |
| D8 | Two-means split replacing the centre | Simple; revisit with R4. |
| D9 | `MAX_CELL_RADIUS = 2.0` | Guards the stereographic scale; lower it for R3. |
| D10 | Sphere only in v0.1 | Matches the public 1992 code; ellipsoid is R5. |
| D11 | Blue Oak Model License 1.0.0 | Permissive; the 1992 data files are excluded and documented. |
