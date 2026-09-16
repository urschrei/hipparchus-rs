use geo::{Distance, Haversine, Point};
use hipparchus_rs::{CellStructure, EARTH_RADIUS, Frame, V3, algebra, build, ddj, hull};

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

#[test]
fn arc_agrees_with_geo_haversine() {
    let mut rng = Rng(7);
    for _ in 0..1000 {
        let (a, b) = (rng.sphere(), rng.sphere());
        let h = Haversine.distance(Point::from(a), Point::from(b));
        assert!((algebra::arc(a, b) * EARTH_RADIUS - h).abs() < 1e-6);
    }
}

type Backend = fn(&[V3]) -> Result<Vec<[u32; 3]>, hull::HullError>;
const BACKENDS: [(&str, Backend); 2] = [
    ("robust", hull::sphere_delaunay),
    ("geo", hull::sphere_delaunay_geo),
];

fn assert_delaunay(pts: &[V3], tris: &[[u32; 3]], name: &str) {
    assert_eq!(tris.len(), 2 * pts.len() - 4, "{name}: Euler");
    // Outward CCW, and empty circumcircle == no site strictly outside a face.
    let c = |v: V3| robust::Coord3D {
        x: v.x,
        y: v.y,
        z: v.z,
    };
    for &[a, b, d] in tris {
        let (a, b, d) = (pts[a as usize], pts[b as usize], pts[d as usize]);
        assert!(
            pts.iter()
                .all(|&p| robust::orient3d(c(a), c(b), c(d), c(p)) >= 0.0),
            "{name}"
        );
    }
}

fn edges(tris: &[[u32; 3]]) -> std::collections::BTreeSet<(u32, u32)> {
    tris.iter()
        .flat_map(|&[a, b, c]| [(a, b), (b, c), (c, a)])
        .collect()
}

#[test]
fn delaunay_is_valid() {
    let mut rng = Rng(42);
    for n in [4, 5, 50, 500, 3000] {
        let pts: Vec<V3> = (0..n).map(|_| rng.sphere()).collect();
        let [(_, a), (_, b)] = BACKENDS.map(|(name, f)| {
            let t = f(&pts).unwrap();
            assert_delaunay(&pts, &t, name);
            (name, t)
        });
        assert_eq!(edges(&a), edges(&b), "backends disagree at n = {n}");
    }
    // Clustered sites (a 1 km patch plus a global scatter) stress the
    // projection: the patch lands far out in the plane when site 0 is in it.
    let mut pts: Vec<V3> = (0..200)
        .map(|_| (V3::from_lat_lon(53.3, -6.3) + rng.sphere() * 1.5e-4).normalized())
        .collect();
    pts.extend((0..200).map(|_| rng.sphere()));
    for (name, f) in BACKENDS {
        assert_delaunay(&pts, &f(&pts).unwrap(), name);
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
    for (name, f) in BACKENDS {
        assert_delaunay(&ico, &f(&ico).unwrap(), name);
        let mut dup = ico.clone();
        dup.push(dup[3]);
        assert_eq!(f(&dup), Err(hull::HullError::Duplicate(12)), "{name}");
    }
}

#[test]
fn walk_finds_nearest_centre() {
    let mut rng = Rng(3);
    let cs = CellStructure::new((0..800).map(|_| rng.sphere()).collect()).unwrap();
    for _ in 0..5000 {
        let p = rng.sphere();
        let brute = (0..cs.len() as u32)
            .max_by(|&a, &b| p.dot(cs.centre(a)).total_cmp(&p.dot(cs.centre(b))))
            .unwrap();
        let start = (rng.f() * cs.len() as f64) as u32;
        assert_eq!(cs.home_cell(p, start).0, brute);
        // and the owning cell's radius really bounds the point
        assert!(algebra::arc(p, cs.centre(brute)) <= cs.radius(brute) * (1.0 + 1e-12));
    }
}

#[test]
fn encode_decode_and_cap_queries() {
    let mut rng = Rng(11);
    let cs = CellStructure::new(build::fibonacci(300)).unwrap();
    let mut hint = 0;
    for _ in 0..5000 {
        let p = rng.sphere();
        let l = cs.encode(p, hint);
        hint = l.cell;
        let err = algebra::arc(p, cs.decode(l));
        assert!(err <= 2.0 * cs.resolution(l.cell), "{err}");
    }
    for _ in 0..50 {
        let (q, r) = (rng.sphere(), rng.f() * 0.5);
        let got = cs.cells_in_cap(q, r, 0);
        for _ in 0..500 {
            let p = rng.sphere();
            if algebra::arc(p, q) <= r {
                assert!(got.contains(&cs.home_cell(p, 0).0));
            }
        }
    }
}
