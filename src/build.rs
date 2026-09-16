//! Cell structures whose cell sizes follow the density of a set of points.
//!
//! [`fibonacci`] gives centres spread almost evenly over the sphere, and so
//! cells of almost equal size.
//!
//! [`adaptive`] starts from Fibonacci centres. It splits each cell that holds
//! more than a given number of points into two, then builds the structure
//! again, until no cell holds too many points. Dense areas get small cells.
//! The seed centres keep the cells in empty areas below
//! [`MAX_CELL_RADIUS`](crate::cells::MAX_CELL_RADIUS).
//!
//! The builder cannot split a cell whose points are identical or less than
//! about 10 cm apart, so such a cell can hold more points than the limit.
//! [`Report::load`] gives the number of points in each cell.
//!
//! Each round builds the whole structure again, so the builder suits data that
//! changes in batches. Section 6 of the [design notes] describes the split.
//!
//! [design notes]: https://github.com/urschrei/hipparchus-rs/blob/main/DESIGN.md#6-adaptive-builder

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
///
/// # Examples
///
/// ```
/// use hipparchus_rs::{CellId, V3, build};
///
/// // 1,600 points on a grid over Ireland and Great Britain.
/// let pts: Vec<V3> = (0..40)
///     .flat_map(|i| (0..40).map(move |j| (50.0 + 0.2 * i as f64, -10.0 + 0.3 * j as f64)))
///     .map(|(lat, lon)| V3::from_lat_lon(lat, lon))
///     .collect();
/// let (cells, report) = build::adaptive(&pts, 100, 12)?;
/// assert!(report.load.iter().all(|&n| n <= 100));
/// assert_eq!(report.load.iter().sum::<usize>(), pts.len());
///
/// // Cells over the points are smaller than cells over the empty Pacific.
/// let (dense, _) = cells.home_cell(pts[0], CellId(0));
/// let (empty, _) = cells.home_cell(V3::from_lat_lon(-40.0, -120.0), dense);
/// assert!(cells.radius(dense) < cells.radius(empty));
/// # Ok::<(), hipparchus_rs::BuildError>(())
/// ```
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
