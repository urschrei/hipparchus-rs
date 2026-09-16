//! Store points in cell-local form and find the points near a place.
//!
//! Builds cells over the 1992 coastline vertices, stores each vertex as a
//! 12-byte `Local` in its cell, then decodes only the cells that the query
//! circle can touch.
//!
//! Usage: cargo run --release --example nearby -- [lat lon radius_km]
//! (default: Dublin, 500 km)

use hipparchus_rs::{CellId, EARTH_RADIUS, Local, V3, arc, build, ddj};
use std::{fs, path::Path};

fn main() {
    let args: Vec<f64> = std::env::args()
        .skip(1)
        .map(|a| a.parse().expect("arguments must be numbers"))
        .collect();
    let [lat, lon, km] = match args[..] {
        [] => [53.35, -6.26, 500.0],
        [lat, lon, km] => [lat, lon, km],
        _ => panic!("usage: nearby [lat lon radius_km]"),
    };

    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let idx = ddj::read_index(&fs::read(dir.join("COAST0.IDX")).unwrap());
    let dat = fs::read(dir.join("COAST.DAT")).unwrap();
    let pts: Vec<V3> = idx
        .iter()
        .flat_map(|rec| ddj::decode_segment(rec, &dat))
        .collect();

    // At most 256 vertices per cell; 120 seed cells cover the empty ocean.
    let (cells, _) = build::adaptive(&pts, 256, 120).expect("coastline cells");
    let mut store: Vec<Vec<Local>> = vec![Vec::new(); cells.len()];
    let mut hint = CellId(0);
    for &p in &pts {
        let local = cells.encode(p, hint);
        store[local.cell.index()].push(local);
        hint = local.cell;
    }
    println!(
        "stored {} vertices in {} cells ({} bytes)",
        pts.len(),
        cells.len(),
        pts.len() * size_of::<Local>()
    );

    let (place, radius) = (V3::from_lat_lon(lat, lon), km * 1e3 / EARTH_RADIUS);
    let candidates = cells.cells_in_cap(place, radius, CellId(0));
    let decoded = candidates
        .iter()
        .map(|&c| store[c.index()].len())
        .sum::<usize>();
    let near: Vec<V3> = candidates
        .iter()
        .flat_map(|&c| &store[c.index()])
        .map(|&local| cells.decode(local))
        .filter(|&p| arc(p, place) <= radius)
        .collect();
    println!(
        "{} vertices within {km} km of ({lat}, {lon}); decoded {decoded} from {} cells",
        near.len(),
        candidates.len()
    );
    if let Some(nearest) = near.iter().map(|&p| arc(p, place)).min_by(f64::total_cmp) {
        println!("nearest vertex: {:.1} km", nearest * EARTH_RADIUS / 1e3);
    }
}
