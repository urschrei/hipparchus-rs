//! Vector algebra on direction cosines, and the stereographic frame that
//! cell-local coordinates use.
//!
//! A position is a unit [`V3`]. [`V3::from_lat_lon`] and [`V3::normalized`]
//! return unit vectors; [`V3::new`] and the arithmetic operators do not
//! normalise.
//!
//! [`arc`] gives the great-circle angle between two positions, with good
//! conditioning at all separations. [`arc_ddj`] reproduces the approximation
//! in the 1992 code. It is never longer than the arc, and is much shorter for
//! distant points (2.33 against π at the antipode), so use it only to exclude
//! candidates.
//!
//! A [`Frame`] projects the sphere stereographically, from the antipode of a
//! centre onto the plane tangent at that centre. A point at arc `r` from the
//! centre maps to plane distance [`plane_extent`]`(r)`, which is at least `r`,
//! so a plane distance is never less than the ground distance it represents.
//!
//! The module is a port of `ALGEBRA.C` from Lukatela and Russell (1992).
//! Unlike the original, [`Frame`] works at the poles. Section 4 of the
//! [design notes] gives the derivations.
//!
//! [design notes]: https://github.com/urschrei/hipparchus-rs/blob/main/DESIGN.md#4-core-representation

use std::ops::{Add, Mul, Neg, Sub};

/// A 3-vector. A unit `V3` is a point on the sphere, given by its direction
/// cosines.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct V3 {
    /// Component towards latitude 0°, longitude 0°.
    pub x: f64,
    /// Component towards latitude 0°, longitude 90° E.
    pub y: f64,
    /// Component towards the north pole.
    pub z: f64,
}

impl V3 {
    /// The vector with these components, not normalised.
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    /// Scalar product.
    pub fn dot(self, o: Self) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }
    /// Vector product.
    pub fn cross(self, o: Self) -> Self {
        Self::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }
    /// Euclidean length.
    pub fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }
    /// This vector scaled to unit length. The zero vector gives NaN.
    pub fn normalized(self) -> Self {
        self * (1.0 / self.norm())
    }
    /// From geodetic-style degrees (treated as geocentric, as in the original).
    pub fn from_lat_lon(lat: f64, lon: f64) -> Self {
        let (lat, lon) = (lat.to_radians(), lon.to_radians());
        Self::new(lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin())
    }
    /// (lat, lon) in degrees.
    pub fn lat_lon(self) -> (f64, f64) {
        (
            self.z.atan2(self.x.hypot(self.y)).to_degrees(),
            self.y.atan2(self.x).to_degrees(),
        )
    }
}

impl Add for V3 {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}
impl Sub for V3 {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}
impl Mul<f64> for V3 {
    type Output = Self;
    fn mul(self, k: f64) -> Self {
        Self::new(self.x * k, self.y * k, self.z * k)
    }
}
impl Neg for V3 {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl From<geo::Point<f64>> for V3 {
    fn from(p: geo::Point<f64>) -> Self {
        Self::from_lat_lon(p.y(), p.x())
    }
}
impl From<V3> for geo::Point<f64> {
    fn from(v: V3) -> Self {
        let (lat, lon) = v.lat_lon();
        geo::Point::new(lon, lat)
    }
}

/// Great-circle angle in radians; well conditioned for all separations.
pub fn arc(a: V3, b: V3) -> f64 {
    a.cross(b).norm().atan2(a.dot(b))
}

/// The DDJ `ArcDist`: `chord + chord³/24`. Underestimates long arcs
/// (2.33 vs π at the antipode), which is safe for culling only.
pub fn arc_ddj(a: V3, b: V3) -> f64 {
    let c2 = (b - a).dot(b - a);
    let c = c2.sqrt();
    c + c2 * c / 24.0
}

/// Stereographic plane extent (unit sphere) of a point `r` radians from the
/// tangency point: `2·tan(r/2)`, which is always ≥ `r`.
pub fn plane_extent(r: f64) -> f64 {
    2.0 * (0.5 * r).tan()
}

/// Tangent-plane frame for the stereographic projection centred on `centre`,
/// projecting from the antipode (as `MapStereo`/`UnMapStereo`).
///
/// The original divides by `hypot(x, y)` and fails at the poles. Here a pole
/// gets a fixed frame (east = +Y). The frame is a deterministic function of
/// the centre, so data encoded against a stored centre always round-trips.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Frame {
    centre: V3,
    east: V3,
    north: V3,
}

impl Frame {
    /// The frame tangent at `centre`, which must be a unit vector.
    pub fn new(centre: V3) -> Self {
        let cp = centre.x.hypot(centre.y);
        let east = if cp > 1e-300 {
            V3::new(-centre.y / cp, centre.x / cp, 0.0)
        } else {
            V3::new(0.0, 1.0, 0.0)
        };
        Self {
            centre,
            east,
            north: centre.cross(east),
        }
    }

    /// The tangent point.
    pub fn centre(&self) -> V3 {
        self.centre
    }

    /// Sphere → plane. Singular only at the antipode of the centre.
    pub fn map(&self, p: V3) -> [f64; 2] {
        let k = 2.0 / (1.0 + self.centre.dot(p));
        [k * self.east.dot(p), k * self.north.dot(p)]
    }

    /// Plane → sphere (closed-form inverse stereographic).
    pub fn unmap(&self, [e, n]: [f64; 2]) -> V3 {
        let r2 = e * e + n * n;
        (self.centre * (4.0 - r2) + (self.east * e + self.north * n) * 4.0).normalized()
    }
}
