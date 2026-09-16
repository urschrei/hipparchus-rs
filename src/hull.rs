//! Delaunay triangulation of points on the sphere.
//!
//! The Delaunay triangulation of unit vectors is their convex hull. Both
//! functions in this module return it as triangles of indices into the input,
//! counter-clockwise seen from outside the sphere:
//!
//! - [`sphere_delaunay_geo`] projects the sites stereographically,
//!   triangulates them in the plane with `geo`, and corrects the result with
//!   exact 3D edge flips. Its time is close to linear in the number of sites.
//!   [`CellStructure`](crate::CellStructure) uses it.
//! - [`sphere_delaunay`] builds the convex hull one site at a time with exact
//!   `orient3d` tests. Its time is O(n²), so use it as a reference for small
//!   inputs.
//!
//! Both functions give the same triangles, except where four or more sites lie
//! exactly on one circle. Sites must be at least [`MIN_SEPARATION`] apart, and
//! both functions identify the same site when they are not. Section 5.2 of the
//! [design notes] gives the details.
//!
//! [design notes]: https://github.com/urschrei/hipparchus-rs/blob/main/DESIGN.md#52-delaunay-triangulation

use crate::algebra::{Frame, V3};
use geo::{MultiPoint, Point, TriangulateDelaunayUnconstrained};
use robust::{Coord3D, orient3d};
use spade::Triangulation;
use spade::handles::FixedVertexHandle;
use std::collections::{HashMap, HashSet};
use std::fmt;

/// Sites must be at least this far apart, in radians (0.64 m on the Earth).
pub const MIN_SEPARATION: f64 = 1e-7;

/// Error from [`sphere_delaunay`] and [`sphere_delaunay_geo`].
///
/// For [`Duplicate`](Self::Duplicate) and [`TooClose`](Self::TooClose), `i` is
/// the lowest index of a site that equals, or is closer than
/// [`MIN_SEPARATION`] to, a site with a lower index. Both functions report
/// the same index.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HullError {
    /// Fewer than four distinct, non-coplanar sites.
    TooFewSites,
    /// Site `i` equals a site with a lower index.
    Duplicate(usize),
    /// Site `i` is closer than [`MIN_SEPARATION`] to a site with a lower index.
    TooClose(usize),
    /// The planar triangulator rejected a projected coordinate: spade accepts
    /// finite magnitudes up to 2^201.
    Projection,
}

impl fmt::Display for HullError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TooFewSites => write!(f, "fewer than four distinct, non-coplanar sites"),
            Self::Duplicate(i) => write!(f, "site {i} equals an earlier site"),
            Self::TooClose(i) => write!(f, "site {i} is too close to an earlier site"),
            Self::Projection => write!(f, "a projected site coordinate is out of range"),
        }
    }
}

impl std::error::Error for HullError {}

fn c3(v: V3) -> Coord3D<f64> {
    Coord3D {
        x: v.x,
        y: v.y,
        z: v.z,
    }
}

/// `true` if `p` is strictly outside the face `(a, b, c)`, which is CCW
/// seen from outside. `orient3d` is positive when `p` is *below* the face.
fn sees(pts: &[V3], [a, b, c]: [u32; 3], p: V3) -> bool {
    orient3d(
        c3(pts[a as usize]),
        c3(pts[b as usize]),
        c3(pts[c as usize]),
        c3(p),
    ) < 0.0
}

/// Triangles (CCW seen from outside) of the Delaunay triangulation of `pts`.
/// Cospherical ties are broken arbitrarily but consistently.
///
/// # Errors
///
/// [`HullError::TooFewSites`], [`HullError::Duplicate`] or
/// [`HullError::TooClose`] as described there.
pub fn sphere_delaunay(pts: &[V3]) -> Result<Vec<[u32; 3]>, HullError> {
    if pts.len() < 4 {
        return Err(HullError::TooFewSites);
    }
    if let Some(e) = separation_error(pts) {
        return Err(e);
    }
    // Three separated sphere points are not collinear, so a fourth point off
    // their plane completes the first tetrahedron.
    let (a, b, c) = (c3(pts[0]), c3(pts[1]), c3(pts[2]));
    let i3 = (3..pts.len())
        .find(|&i| orient3d(a, b, c, c3(pts[i])) != 0.0)
        .ok_or(HullError::TooFewSites)?;

    let t = [0, 1, 2, i3 as u32];
    let mut faces: Vec<[u32; 3]> = [[0, 1, 2], [0, 3, 1], [1, 3, 2], [2, 3, 0]]
        .iter()
        .map(|f| [t[f[0]], t[f[1]], t[f[2]]])
        .collect();
    if sees(pts, faces[0], pts[i3]) {
        // wrong handedness: flip every face
        for f in &mut faces {
            f.swap(1, 2);
        }
    }

    let mut visible_edges = HashSet::new();
    for (i, &p) in pts.iter().enumerate() {
        if i < 3 || i == i3 {
            continue;
        }
        let (vis, keep): (Vec<[u32; 3]>, Vec<[u32; 3]>) =
            faces.iter().partition(|&&f| sees(pts, f, p));
        if vis.is_empty() {
            // Rounding put the site inside the hull. Separated sites do not
            // cause this.
            return Err(HullError::TooClose(i));
        }
        visible_edges.clear();
        for &[a, b, c] in &vis {
            visible_edges.extend([(a, b), (b, c), (c, a)]);
        }
        faces = keep;
        // A visible edge whose twin is not visible is on the horizon; the new
        // face over it keeps the old face's orientation.
        for &(a, b) in &visible_edges {
            if !visible_edges.contains(&(b, a)) {
                faces.push([a, b, i as u32]);
            }
        }
    }
    Ok(faces)
}

/// Same result via `geo`'s Delaunay triangulation: project stereographically
/// from site 0 (circles stay circles), triangulate the other sites in the
/// plane, then join site 0 to the planar convex hull.
///
/// spade's predicates are exact on the projected coordinates. Exact 3D edge
/// flips then correct the result, so it is Delaunay on the sites themselves;
/// only exactly cospherical ties can break differently from
/// [`sphere_delaunay`].
///
/// # Errors
///
/// Any [`HullError`], as described there.
pub fn sphere_delaunay_geo(pts: &[V3]) -> Result<Vec<[u32; 3]>, HullError> {
    if pts.len() < 4 {
        return Err(HullError::TooFewSites);
    }
    // Insert in Morton order: spade walks from the previous insertion, so
    // spatially coherent input keeps each walk short.
    let mut order: Vec<u32> = (1..pts.len() as u32).collect();
    order.sort_by_cached_key(|&i| morton(pts[i as usize]));

    // Tangent at the antipode of site 0 = projection from site 0. The frame
    // is right-handed seen from outside, so planar CCW == spherical CCW.
    let frame = Frame::new(-pts[0]);
    let plane: MultiPoint = order
        .iter()
        .map(|&i| Point::from(frame.map(pts[i as usize]).map(snap_to_zero)))
        .collect();
    // Sites too close to one another cause the failures below; other failures
    // are projection failures.
    let close_or_projection = || separation_error(pts).unwrap_or(HullError::Projection);
    // Sites near site 0 do not project accurately, so check them here.
    if pts[1..].iter().any(|&p| too_close(p, pts[0])) {
        return Err(close_or_projection());
    }
    let dt = plane
        .unconstrained_triangulation_raw()
        .map_err(|_| close_or_projection())?;
    // spade vertex k is the k-th inserted site
    let id = |v: FixedVertexHandle| order[v.index()];
    let short = |a: u32, b: u32| too_close(pts[a as usize], pts[b as usize]);
    // Each site's nearest neighbour is a Delaunay neighbour, so other sites
    // that are too close give a merged vertex or a short planar edge.
    let merged = dt.num_vertices() != order.len();
    if merged
        || dt.undirected_edges().any(|e| {
            let [a, b] = e.vertices().map(|v| id(v.fix()));
            short(a, b)
        })
    {
        return Err(close_or_projection());
    }
    if dt.all_vertices_on_line() {
        return Err(HullError::TooFewSites);
    }
    let mut tris: Vec<[u32; 3]> = dt
        .inner_faces()
        .map(|f| f.vertices().map(|v| id(v.fix())))
        .collect();
    // Hull edges run clockwise around the planar triangulation, so the fan
    // around the projection pole uses them in their given direction.
    tris.extend(
        dt.convex_hull()
            .map(|e| [id(e.from().fix()), id(e.to().fix()), 0]),
    );

    // spade's predicates are exact on the projected coordinates only. Find the
    // edges where a site sees the face across the edge in the exact 3D test.
    // The outer planar face is the fan around site 0.
    let site = |v: Option<FixedVertexHandle>| v.map_or(0, id);
    let mut reflex: Vec<(u32, u32)> = dt
        .undirected_edges()
        .filter_map(|e| {
            let d = e.as_directed();
            let [a, b] = [d.from(), d.to()].map(|v| id(v.fix()));
            let [left, right] = [d, d.rev()].map(|d| site(d.opposite_vertex().map(|v| v.fix())));
            sees(pts, [a, b, left], pts[right as usize]).then_some((a, b))
        })
        .collect();
    // Consecutive hull vertices u, v, w give the fan faces (u, v, 0) and
    // (v, w, 0), which share the edge (v, 0).
    let hull: Vec<u32> = dt.convex_hull().map(|e| id(e.from().fix())).collect();
    for (k, &v) in hull.iter().enumerate() {
        let u = hull[(k + hull.len() - 1) % hull.len()];
        let w = hull[(k + 1) % hull.len()];
        if sees(pts, [u, v, 0], pts[w as usize]) {
            reflex.push((v, 0));
        }
    }
    if !reflex.is_empty() {
        flip_reflex_edges(pts, &mut tris, reflex);
    }
    Ok(tris)
}

/// Flips edges of `tris`, starting from the directed edges in `stack`, until
/// no site sees the face across an edge. Each flip adds a tetrahedron to the
/// enclosed solid, so the flips stop.
fn flip_reflex_edges(pts: &[V3], tris: &mut [[u32; 3]], mut stack: Vec<(u32, u32)>) {
    let mut face = HashMap::with_capacity(3 * tris.len());
    for (t, &[a, b, c]) in tris.iter().enumerate() {
        face.extend([((a, b), t), ((b, c), t), ((c, a), t)]);
    }
    let third = |[x, y, z]: [u32; 3], a: u32, b: u32| {
        if x != a && x != b {
            x
        } else if y != a && y != b {
            y
        } else {
            z
        }
    };
    while let Some((a, b)) = stack.pop() {
        let (Some(&t1), Some(&t2)) = (face.get(&(a, b)), face.get(&(b, a))) else {
            continue;
        };
        // Faces (a, b, c) and (b, a, d) become (a, d, c) and (d, b, c). The
        // flip must not repeat the edge (c, d).
        let (c, d) = (third(tris[t1], a, b), third(tris[t2], a, b));
        if face.contains_key(&(c, d)) || !sees(pts, [a, b, c], pts[d as usize]) {
            continue;
        }
        tris[t1] = [a, d, c];
        tris[t2] = [d, b, c];
        face.remove(&(a, b));
        face.remove(&(b, a));
        face.extend([((a, d), t1), ((d, c), t1), ((b, c), t2), ((c, d), t2)]);
        stack.extend([(a, d), (d, b), (b, c), (c, a)]);
    }
}

/// spade rejects non-zero magnitudes below `MIN_ALLOWED_VALUE` (2^-142). Such
/// a value is smaller than the rounding error of `Frame::map`, so use zero.
fn snap_to_zero(c: f64) -> f64 {
    if c.abs() < spade::MIN_ALLOWED_VALUE {
        0.0
    } else {
        c
    }
}

/// The error for the lowest site index that equals, or is closer than
/// `MIN_SEPARATION` to, a lower index.
fn separation_error(pts: &[V3]) -> Option<HullError> {
    // Sites closer than MIN_SEPARATION are in the same or adjacent boxes.
    let key = |p: V3| [p.x, p.y, p.z].map(|c| (c / MIN_SEPARATION).floor() as i64);
    let mut boxes: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    for (i, &p) in pts.iter().enumerate() {
        let [x, y, z] = key(p);
        let adjacent = (x - 1..=x + 1).flat_map(|x| {
            (y - 1..=y + 1).flat_map(move |y| (z - 1..=z + 1).map(move |z| [x, y, z]))
        });
        let mut close = false;
        for j in adjacent.filter_map(|k| boxes.get(&k)).flatten() {
            if pts[*j] == p {
                return Some(HullError::Duplicate(i));
            }
            close |= too_close(pts[*j], p);
        }
        if close {
            return Some(HullError::TooClose(i));
        }
        boxes.entry([x, y, z]).or_default().push(i);
    }
    None
}

/// `true` if `a` and `b` are closer than `MIN_SEPARATION`. At that distance
/// the chord equals the arc to within rounding.
fn too_close(a: V3, b: V3) -> bool {
    let d = a - b;
    d.dot(d) < MIN_SEPARATION * MIN_SEPARATION
}

/// Interleaved 21-bit quantisations of a unit vector's components.
fn morton(v: V3) -> u64 {
    let q = |c: f64| ((c + 1.0) * 0.5 * ((1 << 21) - 1) as f64) as u64;
    let (x, y, z) = (q(v.x), q(v.y), q(v.z));
    (0..21).fold(0, |m, b| {
        m | ((x >> b) & 1) << (3 * b)
            | ((y >> b) & 1) << (3 * b + 1)
            | ((z >> b) & 1) << (3 * b + 2)
    })
}
