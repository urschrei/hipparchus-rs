//! Data-density-driven cell structures: split any cell holding more than
//! `capacity` points until none does (cf. Geodyssey's population grids).

use crate::algebra::V3;
use crate::cells::{BuildError, CellId, CellStructure};

/// Roughly uniform seed centres (Fibonacci lattice), so even empty regions
/// are covered by cells of bounded size.
pub fn fibonacci(n: usize) -> Vec<V3> {
    let ga = std::f64::consts::PI * (3.0 - 5f64.sqrt());
    (0..n)
        .map(|i| {
            let z = 1.0 - (2.0 * i as f64 + 1.0) / n as f64;
            let r = (1.0 - z * z).sqrt();
            let t = ga * i as f64;
            V3::new(r * t.cos(), r * t.sin(), z)
        })
        .collect()
}

/// Two-means split of a cell's points; `None` if they cannot be separated.
fn split(pts: &[V3]) -> Option<[V3; 2]> {
    let far = |from: V3| {
        pts.iter()
            .copied()
            .min_by(|a, b| a.dot(from).total_cmp(&b.dot(from)))
    };
    let a = far(mean(pts.iter().copied())?)?;
    let b = far(a)?;
    if a == b {
        return None;
    }
    let mut k = [a, b];
    for _ in 0..8 {
        let nearer_second = |p: &&V3| p.dot(k[0]) < p.dot(k[1]);
        k = [
            mean(pts.iter().filter(|p| !nearer_second(p)).copied())?,
            mean(pts.iter().filter(nearer_second).copied())?,
        ];
    }
    Some(k)
}

/// Normalised sum of `pts`; `None` if the sum is zero.
fn mean(pts: impl Iterator<Item = V3>) -> Option<V3> {
    let sum = pts.fold(V3::default(), |a, b| a + b);
    (sum.norm() > 0.0).then(|| sum.normalized())
}

/// Statistics from [`adaptive`].
#[derive(Clone, Debug)]
pub struct Report {
    /// Number of cell structures built.
    pub rounds: usize,
    /// Number of points in each cell of the final structure.
    pub load: Vec<usize>,
}

/// Builds cells from `seeds` Fibonacci-lattice centres, then splits each cell
/// that holds more than `capacity` of `pts` until no cell needs, or can take,
/// a split.
///
/// # Errors
///
/// [`BuildError`] if a structure does not build: for example, if `seeds` is
/// less than four or a split puts two centres too close together.
pub fn adaptive(
    pts: &[V3],
    capacity: usize,
    seeds: usize,
) -> Result<(CellStructure, Report), BuildError> {
    let mut cs = CellStructure::new(fibonacci(seeds))?;
    let mut rounds = 1;
    loop {
        let mut members: Vec<Vec<V3>> = vec![Vec::new(); cs.len()];
        let mut hint = CellId(0);
        for &p in pts {
            // Consecutive vertices are close, so the last cell is a good start.
            hint = cs.home_cell(p, hint).0;
            members[hint.index()].push(p);
        }
        let mut next = Vec::with_capacity(cs.len());
        for (&c, m) in cs.centres().iter().zip(&members) {
            if m.len() > capacity
                && let Some(pair) = split(m)
            {
                next.extend(pair);
            } else {
                next.push(c);
            }
        }
        if next.len() == cs.len() {
            let load = members.iter().map(Vec::len).collect();
            return Ok((cs, Report { rounds, load }));
        }
        cs = CellStructure::new(next)?;
        rounds += 1;
    }
}
