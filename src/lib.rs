#![doc = include_str!("../README.md")]

pub mod algebra;
pub mod build;
pub mod cells;
pub mod ddj;
pub mod hull;

pub use algebra::{Frame, V3, arc};
pub use cells::{CellStructure, Local};

/// Mean Earth radius (m), matching `geo`'s haversine.
pub const EARTH_RADIUS: f64 = 6_371_008.8;
