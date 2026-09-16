//! Spherical Voronoi cells that follow the density of your data, and compact
//! storage of points relative to those cells.
//!
//! This crate reimplements the Hipparchus geopositioning model of Hrvoje
//! Lukatela from the published papers and the C code that accompanied the
//! 1992 *Dr. Dobb's Journal* article. It works on the sphere, and supports
//! points but not lines or regions. It is experimental: expect the API to
//! change.
//!
//! # Concepts
//!
//! - **Positions** are unit vectors ([`V3`]) whose components are the
//!   direction cosines of the surface normal. Distances come from vector
//!   algebra, so there is no singularity at the poles and no seam at the
//!   antimeridian.
//! - **Cells** are the Voronoi cells of a set of centres: each cell holds the
//!   part of the sphere nearer to its centre than to any other centre. A
//!   [`CellStructure`] stores the neighbours, vertices and radius of each
//!   cell, and a [`CellId`] identifies a cell. [`build::adaptive`] chooses
//!   centres so that no cell holds more than a given number of your points.
//! - **Cell-local form** stores a point as a [`Local`]: its cell and two
//!   32-bit offsets from the cell centre, 12 bytes in total.
//!
//! # Example
//!
//! Build cells over a set of points, store the points in cell-local form, and
//! find the points within 400 km of Dublin.
//!
//! ```
//! use hipparchus_rs::{CellId, EARTH_RADIUS, Local, V3, arc, build};
//!
//! // Latitude and longitude, in degrees.
//! let cities = [
//!     (53.35, -6.26), // Dublin
//!     (55.95, -3.19), // Edinburgh
//!     (51.51, -0.13), // London
//!     (48.86, 2.35),  // Paris
//!     (52.52, 13.40), // Berlin
//!     (40.42, -3.70), // Madrid
//! ];
//! let pts: Vec<V3> = cities
//!     .iter()
//!     .map(|&(lat, lon)| V3::from_lat_lon(lat, lon))
//!     .collect();
//!
//! // Split cells until none holds more than two points. The 12 seed centres
//! // limit the size of cells that hold no points.
//! let (cells, _) = build::adaptive(&pts, 2, 12)?;
//!
//! // Group the stored points by cell. Each search for a point's cell starts
//! // from the cell of the previous point.
//! let mut stored: Vec<Vec<Local>> = vec![Vec::new(); cells.len()];
//! let mut hint = CellId(0);
//! for &p in &pts {
//!     let local = cells.encode(p, hint);
//!     hint = local.cell;
//!     stored[local.cell.index()].push(local);
//! }
//!
//! // Decode only the points in cells that the circle can touch, then test
//! // their distance. Arcs are in radians.
//! let dublin = V3::from_lat_lon(53.35, -6.26);
//! let radius = 400e3 / EARTH_RADIUS;
//! let near: Vec<V3> = cells
//!     .cells_in_cap(dublin, radius, hint)
//!     .into_iter()
//!     .flat_map(|c| &stored[c.index()])
//!     .map(|&local| cells.decode(local))
//!     .filter(|&p| arc(dublin, p) <= radius)
//!     .collect();
//! assert_eq!(near.len(), 2); // Dublin and Edinburgh
//! # Ok::<(), hipparchus_rs::BuildError>(())
//! ```
//!
//! # Conventions
//!
//! - Latitude and longitude are in degrees. [`V3::from_lat_lon`] treats
//!   latitude as geocentric. The conversions to and from `geo::Point` use
//!   longitude as `x` and latitude as `y`.
//! - Arcs, radii and resolutions are angles in radians. Multiply by
//!   [`EARTH_RADIUS`] to get metres.
//! - Neighbour rings, Voronoi vertices and triangles are counter-clockwise
//!   seen from outside the sphere.
//! - A [`CellId`] is the index of the cell's centre in the centres given to
//!   [`CellStructure::new`]. Methods that take a `CellId` panic if it is not
//!   less than [`CellStructure::len`].
//!
//! # Limits
//!
//! - Centres must be at least [`hull::MIN_SEPARATION`] apart: 0.64 m on the
//!   Earth.
//! - A cell radius cannot be more than [`cells::MAX_CELL_RADIUS`] (π/2). Use
//!   enough centres to cover areas with no data.
//! - The encoding error is below 1 mm in cells up to about 2,800 km in
//!   radius, and at most 4.5 mm in the widest cells.
//! - You cannot change a [`CellStructure`]. A different set of centres needs
//!   a new structure, and a [`Local`] decodes correctly only against the
//!   structure that encoded it.
//! - Distances are on a sphere of radius [`EARTH_RADIUS`]. They can differ
//!   from distances on the WGS84 ellipsoid by up to about 0.5%.
//!
//! # Further reading
//!
//! - The [README] compares the crate with fixed global grids, gives results
//!   on the Natural Earth coastlines, and lists the papers.
//! - The [design notes] give the derivations, the tests that verify them, and
//!   the roadmap.
//! - The [examples] build cells over the Natural Earth coastlines, query them,
//!   and time the triangulation.
//!
//! [README]: https://github.com/urschrei/hipparchus-rs/blob/main/README.md
//! [design notes]: https://github.com/urschrei/hipparchus-rs/blob/main/DESIGN.md
//! [examples]: https://github.com/urschrei/hipparchus-rs/tree/main/examples

#![warn(missing_docs, clippy::doc_markdown, clippy::missing_errors_doc)]

pub mod algebra;
pub mod build;
pub mod cells;
pub mod hull;

pub use algebra::{Frame, V3, arc};
pub use cells::{BuildError, CellId, CellStructure, Local};
pub use hull::HullError;

/// Mean Earth radius (m), matching `geo`'s haversine.
pub const EARTH_RADIUS: f64 = 6_371_008.8;

/// Compiles and runs the code blocks in the README.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
pub struct ReadmeDoctests;
