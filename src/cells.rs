//! Hipparchus-style Voronoi cell structure (Lukatela 1987).

use crate::algebra::{Frame, V3, arc, plane_extent};
use crate::hull::{HullError, sphere_delaunay_geo};
use std::collections::VecDeque;

/// Largest local coordinate magnitude used (headroom below `i32::MAX`).
const LOCAL_MAX: f64 = 2.0e9;
/// Refuse cells wider than this (radians). Within it, every Voronoi edge lies
/// in the hemisphere about the centre, so the farthest point of a cell is a
/// vertex. Add more centres to cover sparse regions.
pub const MAX_CELL_RADIUS: f64 = std::f64::consts::FRAC_PI_2;

#[derive(Debug, PartialEq)]
pub enum BuildError {
    Hull(HullError),
    CellTooLarge { cell: u32, radius: f64 },
}

impl From<HullError> for BuildError {
    fn from(e: HullError) -> Self {
        Self::Hull(e)
    }
}

/// A point in cell-local form: 4 + 4 + 4 bytes on disk (the paper used a
/// 16-bit cell id; `u32` just raises the cell-count ceiling).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Local {
    pub cell: u32,
    pub e: i32,
    pub n: i32,
}

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
    pub fn new(centres: Vec<V3>) -> Result<Self, BuildError> {
        let tris = sphere_delaunay_geo(&centres)?;
        let n = centres.len();

        // For each site, the (next, prev) pair of each incident triangle,
        // plus that triangle's circumcentre (a Voronoi vertex).
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

    pub fn len(&self) -> usize {
        self.centres.len()
    }
    pub fn is_empty(&self) -> bool {
        self.centres.is_empty()
    }
    pub fn centre(&self, c: u32) -> V3 {
        self.centres[c as usize]
    }
    pub fn radius(&self, c: u32) -> f64 {
        self.radius[c as usize]
    }
    pub fn neighbours(&self, c: u32) -> &[u32] {
        &self.neighbours[c as usize]
    }
    pub fn vertices(&self, c: u32) -> &[V3] {
        &self.vertices[c as usize]
    }
    /// Ground distance represented by one local unit, in radians.
    pub fn resolution(&self, c: u32) -> f64 {
        1.0 / self.scale[c as usize]
    }

    /// Greedy walk from `start` to the cell owning `p`. Returns the cell
    /// and the number of moves. Exact on a Delaunay graph; uses the dot
    /// product as a monotone proxy for distance.
    pub fn home_cell(&self, p: V3, start: u32) -> (u32, usize) {
        let (mut cur, mut steps) = (start, 0);
        let mut best = p.dot(self.centre(cur));
        loop {
            let mut moved = false;
            for &nb in self.neighbours(cur) {
                let d = p.dot(self.centre(nb));
                if d > best {
                    (best, cur, moved) = (d, nb, true);
                }
            }
            if !moved {
                return (cur, steps);
            }
            steps += 1;
        }
    }

    /// Cells that may intersect the cap of `radius` (radians) around `centre`:
    /// a breadth-first search outward from the cap centre's home cell.
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

    pub fn decode(&self, l: Local) -> V3 {
        let s = self.resolution(l.cell);
        self.frames[l.cell as usize].unmap([l.e as f64 * s, l.n as f64 * s])
    }
}
