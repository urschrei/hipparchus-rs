use geo::{Distance, Haversine, Point};
use hegel::TestCase;
use hegel::generators as gs;
use hipparchus_rs::cells::MAX_CELL_RADIUS;
use hipparchus_rs::{CellStructure, EARTH_RADIUS, Frame, V3, algebra, build, ddj, hull};
use std::collections::{BTreeSet, HashMap};

struct Rng(u64);
impl Rng {
    fn f(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn sphere(&mut self) -> V3 {
        let z = 2.0 * self.f() - 1.0;
        let t = std::f64::consts::TAU * self.f();
        let r = (1.0 - z * z).sqrt();
        V3::new(r * t.cos(), r * t.sin(), z)
    }
}

fn rows(csv: &str) -> impl Iterator<Item = Vec<f64>> + '_ {
    csv.lines()
        .map(|l| l.split(',').map(|x| x.parse().unwrap()).collect())
}

fn close(a: V3, b: V3, tol: f64) -> bool {
    (a - b).norm() < tol
}

fn v3([x, y, z]: [f64; 3]) -> V3 {
    V3::new(x, y, z)
}

/// A unit vector, as coordinates. Draws near the origin give the north pole,
/// so failures shrink towards it.
#[hegel::composite]
fn unit_vector(tc: &TestCase) -> [f64; 3] {
    let c = || tc.draw(gs::floats::<f64>().min_value(-1.0).max_value(1.0));
    let v = V3::new(c(), c(), c());
    let u = if v.norm() < 1e-3 {
        V3::new(0.0, 0.0, 1.0)
    } else {
        v.normalized()
    };
    [u.x, u.y, u.z]
}

#[hegel::composite]
fn sites(tc: &TestCase, min: usize, max: usize) -> Vec<[f64; 3]> {
    let n = tc.draw(gs::integers::<usize>().min_value(min).max_value(max));
    tc.draw(gs::vecs(unit_vector()).min_size(n).max_size(n))
}

fn draw_point(tc: &TestCase) -> V3 {
    v3(tc.draw(unit_vector()))
}

/// Sites must be about a metre apart (`HullError::TooClose`).
const MIN_SEPARATION: f64 = 2e-7;

/// Appends each site that is at least `MIN_SEPARATION` from every site
/// already in `out`.
fn separate(mut out: Vec<V3>, sites: impl IntoIterator<Item = V3>) -> Vec<V3> {
    for p in sites {
        if out.iter().all(|&q| algebra::arc(p, q) >= MIN_SEPARATION) {
            out.push(p);
        }
    }
    out
}

fn draw_sites(tc: &TestCase, min: usize, max: usize) -> Vec<V3> {
    separate(Vec::new(), tc.draw(sites(min, max)).into_iter().map(v3))
}

/// A cell structure over a Fibonacci lattice, which keeps cells small, plus
/// arbitrary sites. Site sets that do not build are rejected.
fn draw_structure(tc: &TestCase) -> CellStructure {
    let lattice = tc.draw(gs::integers::<usize>().min_value(8).max_value(200));
    let drawn = tc.draw(sites(0, 200)).into_iter().map(v3);
    CellStructure::new(separate(build::fibonacci(lattice), drawn)).unwrap_or_else(|_| tc.reject())
}

fn draw_cell(tc: &TestCase, cs: &CellStructure) -> u32 {
    tc.draw(gs::integers::<u32>().max_value(cs.len() as u32 - 1))
}

#[test]
fn stereo_matches_original_c() {
    for r in rows(include_str!("data/stereo.csv")) {
        let (p0, pe) = (V3::new(r[0], r[1], r[2]), V3::new(r[3], r[4], r[5]));
        let f = Frame::new(p0);
        let [e, n] = f.map(pe);
        assert!(
            (e - r[6]).abs() < 1e-12 && (n - r[7]).abs() < 1e-12,
            "{r:?}"
        );
        let back = f.unmap([r[6], r[7]]);
        assert!(close(back, V3::new(r[8], r[9], r[10]), 1e-12));
        assert!((algebra::arc_ddj(p0, pe) - r[11]).abs() < 1e-14);
    }
}

#[test]
fn coast_decode_matches_original_c() {
    let idx = ddj::read_index(include_bytes!("data/COAST0.IDX"));
    let dat = include_bytes!("data/COAST.DAT");
    assert_eq!(idx.len(), 992);
    for (rec, r) in idx.iter().zip(rows(include_str!("data/coast0_first.csv"))) {
        assert_eq!(rec.offset as f64, r[0]);
        let first = ddj::decode_segment(rec, dat)[0];
        assert!(close(first, V3::new(r[1], r[2], r[3]), 1e-12));
    }
}

#[test]
fn poles_round_trip() {
    for pole in [V3::new(0.0, 0.0, 1.0), V3::new(0.0, 0.0, -1.0)] {
        let f = Frame::new(pole);
        let p = V3::from_lat_lon(80.0 * pole.z, 123.0);
        assert!(close(f.unmap(f.map(p)), p, 1e-14));
    }
}

#[hegel::test]
fn lat_lon_round_trips(tc: TestCase) {
    let p = draw_point(&tc);
    let (lat, lon) = p.lat_lon();
    assert!(close(V3::from_lat_lon(lat, lon), p, 1e-14));
}

#[hegel::test]
fn arc_agrees_with_geo_haversine(tc: TestCase) {
    let (a, b) = (draw_point(&tc), draw_point(&tc));
    let h = Haversine.distance(Point::from(a), Point::from(b));
    let arc = algebra::arc(a, b);
    // Haversine takes the arcsine of sin(arc / 2), which amplifies rounding
    // by 1 / cos(arc / 2) near the antipode.
    let tolerance = 1e-6 + EARTH_RADIUS * 8.0 * f64::EPSILON / (0.5 * arc).cos();
    assert!((arc * EARTH_RADIUS - h).abs() < tolerance);
}

#[hegel::test]
fn arc_ddj_does_not_exceed_arc(tc: TestCase) {
    let (a, b) = (draw_point(&tc), draw_point(&tc));
    assert!(algebra::arc_ddj(a, b) <= algebra::arc(a, b) + 1e-15);
}

#[hegel::test]
fn frame_round_trips(tc: TestCase) {
    let (c, p) = (draw_point(&tc), draw_point(&tc));
    tc.assume(algebra::arc(c, p) <= MAX_CELL_RADIUS);
    let f = Frame::new(c);
    assert!(close(f.unmap(f.map(p)), p, 1e-14));
}

#[hegel::test]
fn frame_maps_arc_to_plane_extent(tc: TestCase) {
    let (c, p) = (draw_point(&tc), draw_point(&tc));
    let r = algebra::arc(c, p);
    tc.assume(r <= MAX_CELL_RADIUS);
    let [e, n] = Frame::new(c).map(p);
    let want = algebra::plane_extent(r);
    assert!((e.hypot(n) - want).abs() <= 1e-14 * (1.0 + want));
}

/// Largest `orient3d` magnitude that rounding in the projected coordinates
/// can give the wrong sign.
const PROJECTION_ROUNDING: f64 = 1e-14;

type Backend = fn(&[V3]) -> Result<Vec<[u32; 3]>, hull::HullError>;
/// Each backend with the `orient3d` tolerance its predicates give: the geo
/// backend is exact on projected coordinates only.
const BACKENDS: [(&str, Backend, f64); 2] = [
    ("robust", hull::sphere_delaunay, 0.0),
    ("geo", hull::sphere_delaunay_geo, PROJECTION_ROUNDING),
];

fn c3(v: V3) -> robust::Coord3D<f64> {
    robust::Coord3D {
        x: v.x,
        y: v.y,
        z: v.z,
    }
}

fn assert_delaunay(pts: &[V3], tris: &[[u32; 3]], name: &str, tolerance: f64) {
    assert_eq!(tris.len(), 2 * pts.len() - 4, "{name}: Euler");
    // The faces form a closed surface: each directed edge occurs once, and so
    // does its reverse.
    let directed: Vec<(u32, u32)> = tris
        .iter()
        .flat_map(|&[a, b, c]| [(a, b), (b, c), (c, a)])
        .collect();
    let set: BTreeSet<_> = directed.iter().copied().collect();
    assert_eq!(set.len(), directed.len(), "{name}: repeated edge");
    assert!(
        set.iter().all(|&(a, b)| set.contains(&(b, a))),
        "{name}: open edge"
    );
    // Outward CCW, and empty circumcircle == no site strictly outside a face.
    let outside = |o: f64| o < -tolerance;
    for &[a, b, d] in tris {
        let (a, b, d) = (pts[a as usize], pts[b as usize], pts[d as usize]);
        assert!(
            pts.iter()
                .all(|&p| !outside(robust::orient3d(c3(a), c3(b), c3(d), c3(p)))),
            "{name}: site outside a face"
        );
    }
}

fn edges(tris: &[[u32; 3]]) -> BTreeSet<(u32, u32)> {
    tris.iter()
        .flat_map(|&[a, b, c]| [(a, b), (b, c), (c, a)])
        .collect()
}

/// `true` if the site opposite some edge lies on, or within rounding of, the
/// circumcircle of the face across that edge, so that the geo backend can
/// flip the edge.
fn has_near_tie(pts: &[V3], tris: &[[u32; 3]]) -> bool {
    let opposite: HashMap<(u32, u32), u32> = tris
        .iter()
        .flat_map(|&[a, b, c]| [((a, b), c), ((b, c), a), ((c, a), b)])
        .collect();
    let p = |i: u32| c3(pts[i as usize]);
    opposite.iter().any(|(&(a, b), &c)| {
        robust::orient3d(p(a), p(b), p(c), p(opposite[&(b, a)])).abs() <= PROJECTION_ROUNDING
    })
}

#[test]
fn clustered_sites_triangulate() {
    // A 1 km patch plus a global scatter stresses the projection: the patch
    // lands far out in the plane when site 0 is in it.
    let mut rng = Rng(42);
    let mut pts: Vec<V3> = (0..200)
        .map(|_| (V3::from_lat_lon(53.3, -6.3) + rng.sphere() * 1.5e-4).normalized())
        .collect();
    pts.extend((0..200).map(|_| rng.sphere()));
    for (name, f, tolerance) in BACKENDS {
        assert_delaunay(&pts, &f(&pts).unwrap(), name, tolerance);
    }
}

#[hegel::test]
fn delaunay_backends_are_valid(tc: TestCase) {
    let pts = draw_sites(&tc, 0, 300);
    for (name, f, tolerance) in BACKENDS {
        if let Ok(tris) = f(&pts) {
            assert_delaunay(&pts, &tris, name, tolerance);
        }
    }
}

#[hegel::test]
fn delaunay_backends_agree_without_near_ties(tc: TestCase) {
    let pts = draw_sites(&tc, 4, 300);
    let (Ok(exact), Ok(geo)) = (hull::sphere_delaunay(&pts), hull::sphere_delaunay_geo(&pts))
    else {
        tc.reject()
    };
    tc.assume(!has_near_tie(&pts, &exact));
    assert_eq!(edges(&exact), edges(&geo));
}

#[hegel::test]
fn backends_report_the_first_duplicate(tc: TestCase) {
    let mut pts = draw_sites(&tc, 4, 100);
    tc.assume(BACKENDS.iter().all(|(_, f, _)| f(&pts).is_ok()));
    let copies = tc.draw(gs::integers::<usize>().min_value(1).max_value(4));
    for _ in 0..copies {
        let from = tc.draw(gs::integers::<usize>().max_value(pts.len() - 1));
        let to = tc.draw(gs::integers::<usize>().max_value(pts.len()));
        pts.insert(to, pts[from]);
    }
    let first = (1..pts.len()).find(|&i| pts[..i].contains(&pts[i]));
    for (name, f, _) in BACKENDS {
        assert_eq!(
            f(&pts),
            Err(hull::HullError::Duplicate(first.unwrap())),
            "{name}"
        );
    }
}

#[test]
fn cospherical_and_duplicate_sites() {
    // Icosahedron: each vertex's five neighbours are coplanar.
    let g = (1.0 + 5f64.sqrt()) / 2.0;
    let mut ico = Vec::new();
    for s1 in [-1.0, 1.0] {
        for s2 in [-g, g] {
            ico.extend([
                V3::new(0.0, s1, s2),
                V3::new(s1, s2, 0.0),
                V3::new(s2, 0.0, s1),
            ]);
        }
    }
    let ico: Vec<V3> = ico.into_iter().map(V3::normalized).collect();
    let cs = CellStructure::new(ico.clone()).unwrap();
    assert!((0..12).all(|c| cs.neighbours(c).len() == 5));
    for (name, f, tolerance) in BACKENDS {
        assert_delaunay(&ico, &f(&ico).unwrap(), name, tolerance);
        let mut dup = ico.clone();
        dup.push(dup[3]);
        assert_eq!(f(&dup), Err(hull::HullError::Duplicate(12)), "{name}");
    }
}

#[hegel::test]
fn neighbour_rings_are_symmetric(tc: TestCase) {
    let cs = draw_structure(&tc);
    for c in 0..cs.len() as u32 {
        let ring = cs.neighbours(c);
        let distinct: BTreeSet<_> = ring.iter().collect();
        assert_eq!(distinct.len(), ring.len(), "cell {c}: repeated neighbour");
        for &nb in ring {
            assert!(nb != c && cs.neighbours(nb).contains(&c), "{c} and {nb}");
        }
    }
}

#[hegel::test]
fn voronoi_vertices_are_nearest_to_their_sites(tc: TestCase) {
    let cs = draw_structure(&tc);
    for c in 0..cs.len() as u32 {
        let (ring, verts) = (cs.neighbours(c), cs.vertices(c));
        assert_eq!(ring.len(), verts.len());
        for (k, &v) in verts.iter().enumerate() {
            // Vertex k is equidistant from the centre and neighbours k and
            // k + 1, and no centre is nearer.
            let d = v.dot(cs.centre(c));
            for nb in [ring[k], ring[(k + 1) % ring.len()]] {
                assert!((v.dot(cs.centre(nb)) - d).abs() < 1e-12, "{c}, {k}");
            }
            for other in 0..cs.len() as u32 {
                assert!(v.dot(cs.centre(other)) <= d + 1e-12, "{c}, {k}");
            }
        }
    }
}

#[hegel::test]
fn walk_finds_nearest_centre(tc: TestCase) {
    let cs = draw_structure(&tc);
    let (p, start) = (draw_point(&tc), draw_cell(&tc, &cs));
    let nearest = (0..cs.len() as u32)
        .map(|c| p.dot(cs.centre(c)))
        .fold(f64::NEG_INFINITY, f64::max);
    let (home, _) = cs.home_cell(p, start);
    assert_eq!(p.dot(cs.centre(home)), nearest);
}

#[hegel::test]
fn cell_radius_bounds_its_points(tc: TestCase) {
    let cs = draw_structure(&tc);
    let (p, start) = (draw_point(&tc), draw_cell(&tc, &cs));
    let (home, _) = cs.home_cell(p, start);
    assert!(algebra::arc(p, cs.centre(home)) <= cs.radius(home) * (1.0 + 1e-12));
}

#[hegel::test]
fn encoding_error_is_bounded(tc: TestCase) {
    let cs = draw_structure(&tc);
    let (p, hint) = (draw_point(&tc), draw_cell(&tc, &cs));
    let l = cs.encode(p, hint);
    // Rounding moves each plane coordinate by at most half a unit, and one
    // plane unit is at most one unit of ground distance.
    let bound = std::f64::consts::FRAC_1_SQRT_2 * cs.resolution(l.cell);
    let err = algebra::arc(p, cs.decode(l));
    assert!(err <= bound + 1e-15, "{err} > {bound}");
}

#[hegel::test]
fn encoding_does_not_saturate(tc: TestCase) {
    let cs = draw_structure(&tc);
    let (p, hint) = (draw_point(&tc), draw_cell(&tc, &cs));
    let l = cs.encode(p, hint);
    assert!(l.e.unsigned_abs() < i32::MAX as u32 && l.n.unsigned_abs() < i32::MAX as u32);
}

#[hegel::test]
fn cap_query_contains_every_home_cell(tc: TestCase) {
    let cs = draw_structure(&tc);
    let q = draw_point(&tc);
    let r = tc.draw(
        gs::floats::<f64>()
            .min_value(0.0)
            .max_value(std::f64::consts::PI),
    );
    // A point at arc `t <= r` from `q`.
    let t = r * tc.draw(gs::floats::<f64>().min_value(0.0).max_value(1.0));
    let theta = tc.draw(
        gs::floats::<f64>()
            .min_value(0.0)
            .max_value(std::f64::consts::TAU),
    );
    let rho = algebra::plane_extent(t);
    let p = Frame::new(q).unmap([rho * theta.cos(), rho * theta.sin()]);
    tc.assume(algebra::arc(p, q) <= r);
    let got = cs.cells_in_cap(q, r, draw_cell(&tc, &cs));
    let (home, _) = cs.home_cell(p, draw_cell(&tc, &cs));
    assert!(got.contains(&home));
}

#[hegel::test(test_cases = 50)]
fn adaptive_cells_respect_capacity(tc: TestCase) {
    let pts = draw_sites(&tc, 0, 500);
    let capacity = tc.draw(gs::integers::<usize>().max_value(64));
    let seeds = tc.draw(gs::integers::<usize>().min_value(12).max_value(100));
    let (cs, report) = build::adaptive(&pts, capacity, seeds).unwrap_or_else(|_| tc.reject());
    assert_eq!(report.load.len(), cs.len());
    assert_eq!(report.load.iter().sum::<usize>(), pts.len());
    // Only identical points cannot be split, so a cell exceeds its capacity
    // only if it holds a single point.
    assert!(
        report.load.iter().all(|&n| n <= capacity.max(1)),
        "{:?}",
        report.load
    );
}
