//! Reads the coastline lines that the examples use.

use geo::{Geometry, GeometryCollection};
use geojson::GeoJson;
use hipparchus_rs::V3;
use std::{
    fs,
    path::{Path, PathBuf},
};

/// The Natural Earth 1:50m coastlines in `tests/data`.
pub fn default_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/ne_50m_coastline.geojson")
}

/// Reads each line of a GeoJSON file as a list of vertices.
///
/// # Panics
///
/// If the file cannot be read or parsed, or holds a geometry that is not a
/// line.
pub fn load(path: &Path) -> Vec<Vec<V3>> {
    let json = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let geojson: GeoJson = json.parse().unwrap();
    GeometryCollection::<f64>::try_from(&geojson)
        .unwrap()
        .into_iter()
        .flat_map(|g| match g {
            Geometry::LineString(line) => vec![line],
            Geometry::MultiLineString(lines) => lines.0,
            other => panic!("not a line: {other:?}"),
        })
        .map(|line| line.points().map(V3::from).collect())
        .collect()
}
