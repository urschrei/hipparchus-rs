//! Hipparchus-style Voronoi cell structure (Lukatela 1987).

use crate::algebra::{Frame, V3, arc, plane_extent};
use crate::hull::{HullError, sphere_delaunay_geo};
use std::collections::VecDeque;
use std::{error, fmt};

/// Largest local coordinate magnitude used (headroom below `i32::MAX`).
const LOCAL_MAX: f64 = 2.0e9;
/// Refuse cells wider than this (radians). Within it, every Voronoi edge lies
/// in the hemisphere about the centre, so the farthest point of a cell is a
/// vertex. Add more centres to cover sparse regions.
pub const MAX_CELL_RADIUS: f64 = std::f64::consts::FRAC_PI_2;

/// Error from [`CellStructure::new`].
#[derive(Clone, Debug, PartialEq)]
pub enum BuildError {
    /// The centres cannot be triangulated.
    Hull(HullError),
    /// A cell radius exceeds [`MAX_CELL_RADIUS`].
    CellTooLarge {
        /// The cell identifier.
        cell: u32,
        /// The cell radius in radians.
        radius: f64,
    },
}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Hull(_) => write!(f, "cannot triangulate the centres"),
            Self::CellTooLarge { cell, radius } => write!(
                f,
                "cell {cell} has radius {radius} rad, more than MAX_CELL_RADIUS"
            ),
        }
    }
}

impl error::Error for BuildError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Hull(e) => Some(e),
            Self::CellTooLarge { .. } => None,
        }
    }
}

impl From<HullError> for BuildError {
    fn from(e: HullError) -> Self {
        Self::Hull(e)
    }
}

/// A point in cell-local form: a cell identifier and stereographic offsets
/// from the cell centre, 12 bytes in total. AC8 used a 16-bit cell
/// identifier.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Local {
    /// Home cell identifier.
    pub cell: u32,
    /// East offset in local units.
    pub e: i32,
    /// North offset in local units.
    pub n: i32,
}

/// The Voronoi cells of a set of centres, with the neighbour rings, vertices,
/// radii and stereographic frames that queries and point encoding use.
///
/// Methods that take a cell identifier panic if it is not less than
/// [`len`](Self::len).
#[derive(Clone, Debug)]
pub struct CellStructure {
    centres: Vec<V3>,
    /// Delaunay neighbours, circular and CCW seen from outside.
    neighbours: Vec<Vec<u32>>,
    /// Voronoi vertices of each cell, in the same circular order.
    vertices: Vec<Vec<V3>>,
    /// Maximum arc from the centre to any point of the cell.
    radius: Vec<f64>,
    frames: Vec<Frame>,
    /// Local units per unit of stereographic plane distance.
    scale: Vec<f64>,
}

impl CellStructure {
    /// Builds the cells of `centres`, which must be unit vectors about a
    /// metre or more apart. Cell `i` has centre `centres[i]`.
    ///
    /// # Errors
    ///
    /// [`BuildError::Hull`] if the centres cannot be triangulated, and
    /// [`BuildError::CellTooLarge`] if a cell is too wide.
    pub fn new(centres: Vec<V3>) -> Result<Self, BuildError> {
        let tris = sphere_delaunay_geo(&centres)?;
        let n = centres.len();

        // For each site `a`, one `(b, c, circumcentre)` per incident
        // triangle `(a, b, c)`. The circumcentre is a Voronoi vertex.
        let mut fans: Vec<Vec<(u32, u32, V3)>> = vec![Vec::new(); n];
        for &[a, b, c] in &tris {
            let (pa, pb, pc) = (
                centres[a as usize],
                centres[b as usize],
                centres[c as usize],
            );
            let cc = (pb - pa).cross(pc - pa).normalized();
            fans[a as usize].push((b, c, cc));
            fans[b as usize].push((c, a, cc));
            fans[c as usize].push((a, b, cc));
        }

        let mut neighbours = Vec::with_capacity(n);
        let mut vertices = Vec::with_capacity(n);
        let mut radius = Vec::with_capacity(n);
        for (i, fan) in fans.iter().enumerate() {
            // Chain the fan: triangle (i, x, y) is followed by (i, y, z).
            let mut ring = Vec::with_capacity(fan.len());
            let mut verts = Vec::with_capacity(fan.len());
            let mut t = fan[0];
            for _ in 0..fan.len() {
                ring.push(t.0);
                verts.push(t.2);
                t = *fan.iter().find(|u| u.0 == t.1).expect("closed fan");
            }
            let r = verts
                .iter()
                .map(|&v| arc(centres[i], v))
                .fold(0.0, f64::max);
            if r > MAX_CELL_RADIUS {
                return Err(BuildError::CellTooLarge {
                    cell: i as u32,
                    radius: r,
                });
            }
            neighbours.push(ring);
            vertices.push(verts);
            radius.push(r);
        }

        let frames = centres.iter().map(|&c| Frame::new(c)).collect();
        // Fixes the DDJ overflow: size the scale by 2·tan(r/2), not r.
        let scale = radius
            .iter()
            .map(|&r| LOCAL_MAX / plane_extent(r))
            .collect();
        Ok(Self {
            centres,
            neighbours,
            vertices,
            radius,
            frames,
            scale,
        })
    }

    /// Number of cells.
    pub fn len(&self) -> usize {
        self.centres.len()
    }
    /// `true` if there are no cells.
    pub fn is_empty(&self) -> bool {
        self.centres.is_empty()
    }
    /// Identifiers of all cells.
    pub fn ids(&self) -> impl ExactSizeIterator<Item = u32> + use<> {
        0..self.len() as u32
    }
    /// Centres of all cells, indexed by cell identifier.
    pub fn centres(&self) -> &[V3] {
        &self.centres
    }
    /// Centre of cell `c`.
    pub fn centre(&self, c: u32) -> V3 {
        self.centres[c as usize]
    }
    /// Largest arc from the centre of cell `c` to a point of the cell, in
    /// radians.
    pub fn radius(&self, c: u32) -> f64 {
        self.radius[c as usize]
    }
    /// Cells that share an edge with cell `c`, counter-clockwise seen from
    /// outside.
    pub fn neighbours(&self, c: u32) -> &[u32] {
        &self.neighbours[c as usize]
    }
    /// Voronoi vertices of cell `c`. Vertex `k` is between neighbours `k` and
    /// `k + 1`.
    pub fn vertices(&self, c: u32) -> &[V3] {
        &self.vertices[c as usize]
    }
    /// Plane distance of one local unit in cell `c`, in radians. The ground
    /// distance is not more than this.
    pub fn resolution(&self, c: u32) -> f64 {
        1.0 / self.scale[c as usize]
    }

    /// Greedy walk from cell `start` to the cell whose centre is nearest to
    /// `p`. Returns the cell and the number of moves. Correct on a Delaunay
    /// graph; uses the dot product as a monotone proxy for distance.
    pub fn home_cell(&self, p: V3, start: u32) -> (u32, usize) {
        let (mut cur, mut steps) = (start, 0);
        let mut best = p.dot(self.centre(cur));
        loop {
            let prev = cur;
            for &nb in self.neighbours(prev) {
                let d = p.dot(self.centre(nb));
                if d > best {
                    (best, cur) = (d, nb);
                }
            }
            if cur == prev {
                return (cur, steps);
            }
            steps += 1;
        }
    }

    /// Cells that may intersect the cap of `radius` (radians) around `centre`:
    /// a breadth-first search outward from the cap centre's home cell, which
    /// the walk finds from cell `hint`.
    pub fn cells_in_cap(&self, centre: V3, radius: f64, hint: u32) -> Vec<u32> {
        let (home, _) = self.home_cell(centre, hint);
        let mut seen = vec![false; self.len()];
        let mut out = Vec::new();
        let mut queue = VecDeque::from([home]);
        seen[home as usize] = true;
        while let Some(c) = queue.pop_front() {
            if arc(centre, self.centre(c)) > radius + self.radius(c) {
                continue;
            }
            out.push(c);
            for &nb in self.neighbours(c) {
                if !std::mem::replace(&mut seen[nb as usize], true) {
                    queue.push_back(nb);
                }
            }
        }
        out
    }

    /// Encodes `p` in its home cell, walking from cell `hint`. For sequential
    /// points, the cell of the previous point is a good hint.
    pub fn encode(&self, p: V3, hint: u32) -> Local {
        let (cell, _) = self.home_cell(p, hint);
        let [e, n] = self.frames[cell as usize].map(p);
        let s = self.scale[cell as usize];
        Local {
            cell,
            e: (e * s).round() as i32,
            n: (n * s).round() as i32,
        }
    }

    /// The point that `l` encodes.
    pub fn decode(&self, l: Local) -> V3 {
        let s = self.resolution(l.cell);
        self.frames[l.cell as usize].unmap([l.e as f64 * s, l.n as f64 * s])
    }
}
