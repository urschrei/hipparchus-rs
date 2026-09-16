//! Build a density-driven cell structure over the DDJ 1992 world coastlines.
//! Usage: cargo run --release --example coast -- <dir with COAST0.IDX, COAST.DAT, COAST0.LNS>

use hipparchus_rs::{CellId, CellStructure, EARTH_RADIUS, V3, arc, build, ddj};
use std::{fmt::Write, fs, path::Path, time::Instant};

fn main() {
    let dir = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "tests/data".into());
    let dir = Path::new(&dir);
    let idx = ddj::read_index(&fs::read(dir.join("COAST0.IDX")).unwrap());
    let dat = fs::read(dir.join("COAST.DAT")).unwrap();
    let segs: Vec<Vec<V3>> = idx.iter().map(|r| ddj::decode_segment(r, &dat)).collect();
    let pts: Vec<V3> = segs.iter().flatten().copied().collect();
    let km = |r: f64| r * EARTH_RADIUS / 1e3;

    // How lossy was the 1992 16-bit scheme? Compare with the ASCII source.
    if let Ok(lns) = fs::read_to_string(dir.join("COAST0.LNS")) {
        let src: Vec<V3> = lns
            .lines()
            .filter(|l| !l.starts_with('*'))
            .map(|l| {
                let mut it = l.split([' ', ',', '\t']).filter(|s| !s.is_empty());
                let mut f = || it.next().unwrap().parse::<f64>().unwrap();
                V3::from_lat_lon(f(), f())
            })
            .collect();
        assert_eq!(src.len(), pts.len());
        let worst = src
            .iter()
            .zip(&pts)
            .map(|(&a, &b)| arc(a, b))
            .fold(0.0, f64::max);
        println!(
            "DDJ i16 encoding: {} vertices, worst error {:.1} m",
            pts.len(),
            km(worst) * 1e3
        );
    }

    let t = Instant::now();
    let (cs, rep) = build::adaptive(&pts, 256, 120).unwrap();
    let mut radii: Vec<f64> = cs.ids().map(|c| cs.radius(c)).collect();
    radii.sort_by(f64::total_cmp);
    println!(
        "{} cells in {} rounds ({:.2?}); max load {}; cell radius {:.0} / {:.0} / {:.0} km (min/median/max)",
        cs.len(),
        rep.rounds,
        t.elapsed(),
        rep.load.iter().max().unwrap(),
        km(radii[0]),
        km(radii[radii.len() / 2]),
        km(radii[radii.len() - 1]),
    );

    // Encode everything, walking from the previous vertex's cell.
    let t = Instant::now();
    let (mut hint, mut steps, mut worst) = (CellId(0), 0, 0f64);
    for &p in &pts {
        let (c, s) = cs.home_cell(p, hint);
        steps += s;
        let l = cs.encode(p, c);
        worst = worst.max(arc(p, cs.decode(l)));
        hint = l.cell;
    }
    println!(
        "encode+decode {} pts in {:.2?}; {:.3} walk steps/pt; worst error {:.2} mm; 12 bytes/pt",
        pts.len(),
        t.elapsed(),
        steps as f64 / pts.len() as f64,
        km(worst) * 1e6
    );
    let cold: usize = pts
        .iter()
        .step_by(97)
        .map(|&p| cs.home_cell(p, CellId(0)).1)
        .sum();
    println!(
        "cold walks from cell 0: {:.1} steps/pt",
        cold as f64 / pts.len().div_ceil(97) as f64
    );

    // Cap query around Dublin, checked against a full scan.
    let (q, r) = (V3::from_lat_lon(53.35, -6.26), 500e3 / EARTH_RADIUS);
    let cells = cs.cells_in_cap(q, r, CellId(0));
    let mut owner = vec![false; cs.len()];
    cells.iter().for_each(|&c| owner[c.index()] = true);
    let (mut touched, mut found, mut truth, mut h) = (0, 0, 0, CellId(0));
    for &p in &pts {
        h = cs.home_cell(p, h).0;
        let inside = arc(p, q) <= r;
        truth += inside as usize;
        if owner[h.index()] {
            touched += 1;
            found += inside as usize;
        }
    }
    println!(
        "Dublin 500 km cap: {} cells, {} of {} vertices read ({:.1}%), found {}/{}",
        cells.len(),
        touched,
        pts.len(),
        100.0 * touched as f64 / pts.len() as f64,
        found,
        truth
    );

    fs::write("cells.svg", svg(&cs, &segs, V3::from_lat_lon(35.0, 10.0))).unwrap();
    println!("wrote cells.svg");
}

/// Orthographic view: cell boundaries (subdivided great-circle arcs) and coasts.
fn svg(cs: &CellStructure, segs: &[Vec<V3>], view: V3) -> String {
    let (size, rad) = (800.0, 390.0);
    let east = V3::new(0.0, 0.0, 1.0).cross(view).normalized();
    let north = view.cross(east);
    let project = |p: V3| {
        (
            size / 2.0 + rad * p.dot(east),
            size / 2.0 - rad * p.dot(north),
        )
    };
    let mut out = format!(
        "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {size} {size}'>\
         <circle cx='{c}' cy='{c}' r='{rad}' fill='#f4f1ea' stroke='#999'/>",
        c = size / 2.0
    );
    for c in cs.ids() {
        let vs = cs.vertices(c);
        let ring = (0..vs.len())
            .flat_map(|i| {
                let (a, b) = (vs[i], vs[(i + 1) % vs.len()]);
                (0..8).map(move |k| (a * (8 - k) as f64 + b * k as f64).normalized())
            })
            .chain(std::iter::once(vs[0]));
        let style = "stroke='#c0392b' stroke-width='0.4'";
        path(&mut out, ring, view, &project, style);
    }
    for s in segs {
        let style = "stroke='#2c3e50' stroke-width='0.5'";
        path(&mut out, s.iter().copied(), view, &project, style);
    }
    out + "</svg>"
}

/// Appends a path through `pts`, lifting the pen on the far hemisphere.
fn path(
    out: &mut String,
    pts: impl IntoIterator<Item = V3>,
    view: V3,
    project: &impl Fn(V3) -> (f64, f64),
    style: &str,
) {
    let mut d = String::new();
    let mut pen = false;
    for p in pts {
        if p.dot(view) < 0.0 {
            pen = false;
            continue;
        }
        let (x, y) = project(p);
        write!(d, "{}{x:.1},{y:.1}", if pen { 'L' } else { 'M' }).unwrap();
        pen = true;
    }
    write!(out, "<path d='{d}' fill='none' {style}/>").unwrap();
}
