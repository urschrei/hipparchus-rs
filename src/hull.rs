//! Spherical Delaunay triangulation.
//!
//! - [`sphere_delaunay_geo`]: stereographic projection + `geo`'s (spade)
//!   planar Delaunay. Near-linear; used by `CellStructure`.
//! - [`sphere_delaunay`]: incremental 3D convex hull with exact
//!   `orient3d`. O(n²); kept as the reference implementation.

use crate::algebra::V3;
use robust::{Coord3D, orient3d};
use std::collections::{HashMap, HashSet};

#[derive(Debug, PartialEq)]
pub enum HullError {
    /// Fewer than four distinct, non-coplanar sites.
    TooFewSites,
    /// Site `i` is the first site equal to an earlier site.
    Duplicate(usize),
    /// Site `i` fell inside the hull through rounding: it is within about
    /// a metre of another site. Sites must be separated more than that.
    TooClose(usize),
    /// The planar triangulator rejected a projected coordinate: spade accepts
    /// finite magnitudes up to 2^201.
    Projection,
}

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
pub fn sphere_delaunay(pts: &[V3]) -> Result<Vec<[u32; 3]>, HullError> {
    let n = pts.len();
    // Three distinct sphere points are never collinear, so we only need a
    // fourth point off their plane.
    let i1 = (1..n)
        .find(|&i| pts[i] != pts[0])
        .ok_or(HullError::TooFewSites)?;
    let i2 = (i1 + 1..n)
        .find(|&i| pts[i] != pts[0] && pts[i] != pts[i1])
        .ok_or(HullError::TooFewSites)?;
    let (a, b, c) = (c3(pts[0]), c3(pts[i1]), c3(pts[i2]));
    let i3 = (i2 + 1..n)
        .find(|&i| orient3d(a, b, c, c3(pts[i])) != 0.0)
        .ok_or(HullError::TooFewSites)?;

    let t = [0u32, i1 as u32, i2 as u32, i3 as u32];
    let mut faces: Vec<[u32; 3]> = [[0, 1, 2], [0, 3, 1], [1, 3, 2], [2, 3, 0]]
        .iter()
        .map(|f| [t[f[0]], t[f[1]], t[f[2]]])
        .collect();
    if sees(pts, faces[0], pts[i3]) {
        // wrong handedness: flip every face
        faces.iter_mut().for_each(|f| f.swap(1, 2));
    }

    let mut visible_edges = HashSet::new();
    for (i, &p) in pts.iter().enumerate() {
        if i == 0 || i == i1 || i == i2 || i == i3 {
            continue;
        }
        let (vis, keep): (Vec<[u32; 3]>, Vec<[u32; 3]>) =
            faces.iter().partition(|&&f| sees(pts, f, p));
        if vis.is_empty() {
            let dup = pts[..i].contains(&p);
            return Err(if dup {
                HullError::Duplicate(i)
            } else {
                HullError::TooClose(i)
            });
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
/// Predicates are exact on the *projected* coordinates, so cospherical ties
/// may break differently from [`sphere_delaunay`], and sites within rounding
/// of one another are reported as [`HullError::TooClose`].
pub fn sphere_delaunay_geo(pts: &[V3]) -> Result<Vec<[u32; 3]>, HullError> {
    use crate::algebra::Frame;
    use geo::{MultiPoint, Point, TriangulateDelaunayUnconstrained};
    use spade::Triangulation;

    if pts.len() < 4 {
        return Err(HullError::TooFewSites);
    }
    if let Some(i) = first_duplicate(pts) {
        return Err(HullError::Duplicate(i));
    }
    // Insert in Morton order: spade walks from the previous insertion, so
    // spatially coherent input keeps each walk short.
    let mut order: Vec<u32> = (1..pts.len() as u32).collect();
    order.sort_by_cached_key(|&i| morton(pts[i as usize]));

    // Tangent at the antipode of site 0 = projection from site 0. The frame
    // is right-handed seen from outside, so planar CCW == spherical CCW.
    let frame = Frame::new(-pts[0]);
    let mut plane = Vec::with_capacity(order.len());
    for &i in &order {
        let [x, y] = frame.map(pts[i as usize]).map(snap_to_zero);
        plane.push(Point::new(x, y));
    }
    let dt = MultiPoint::new(plane.clone())
        .unconstrained_triangulation_raw()
        .map_err(|_| HullError::Projection)?;
    if dt.num_vertices() != order.len() {
        // spade merged distinct sites with equal projected coordinates:
        // report the lowest index that coincides with an earlier site
        let mut seen = HashMap::new();
        let mut first = u32::MAX;
        for (&i, q) in order.iter().zip(&plane) {
            if let Some(&j) = seen.get(&bits([q.x(), q.y()])) {
                first = first.min(i.max(j));
            } else {
                seen.insert(bits([q.x(), q.y()]), i);
            }
        }
        return Err(HullError::TooClose(first as usize));
    }
    if dt.all_vertices_on_line() {
        return Err(HullError::TooFewSites);
    }
    // spade vertex k is the k-th inserted site
    let id = |v: spade::handles::FixedVertexHandle| order[v.index()];
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
    Ok(tris)
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

/// Bit patterns of `coords`, with `-0.0` as `0.0` so that equal keys mean `==`.
fn bits<const N: usize>(coords: [f64; N]) -> [u64; N] {
    coords.map(|c| (c + 0.0).to_bits())
}

/// Index of the first site equal to an earlier site.
fn first_duplicate(pts: &[V3]) -> Option<usize> {
    let mut seen = HashSet::with_capacity(pts.len());
    pts.iter().position(|p| !seen.insert(bits([p.x, p.y, p.z])))
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
