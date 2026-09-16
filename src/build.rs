//! Data-density-driven cell structures: split any cell holding more than
//! `capacity` points until none does (cf. Geodyssey's population grids).

use crate::algebra::V3;
use crate::cells::{BuildError, CellStructure};

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
    let mean = |v: &mut dyn Iterator<Item = V3>| {
        let s = v.fold(V3::default(), |a, b| a + b);
        (s.norm() > 0.0).then(|| s.normalized())
    };
    let c = mean(&mut pts.iter().copied())?;
    let far = |from: V3| {
        *pts.iter()
            .min_by(|a, b| a.dot(from).total_cmp(&b.dot(from)))
            .unwrap()
    };
    let mut k = [far(c), V3::default()];
    k[1] = far(k[0]);
    if k[0] == k[1] {
        return None;
    }
    for _ in 0..8 {
        let side = |p: &V3| p.dot(k[0]) < p.dot(k[1]);
        k = [
            mean(&mut pts.iter().filter(|p| !side(p)).copied())?,
            mean(&mut pts.iter().filter(|p| side(p)).copied())?,
        ];
    }
    Some(k)
}

#[derive(Clone, Debug)]
pub struct Report {
    pub rounds: usize,
    pub load: Vec<usize>,
}

pub fn adaptive(
    pts: &[V3],
    capacity: usize,
    seeds: usize,
) -> Result<(CellStructure, Report), BuildError> {
    let mut centres = fibonacci(seeds);
    for rounds in 1.. {
        let cs = CellStructure::new(centres.clone())?;
        let mut members: Vec<Vec<V3>> = vec![Vec::new(); cs.len()];
        let mut hint = 0;
        for &p in pts {
            // Consecutive vertices are close, so the last cell is a good start.
            hint = cs.home_cell(p, hint).0;
            members[hint as usize].push(p);
        }
        let mut split_any = false;
        let mut next = Vec::with_capacity(centres.len());
        for (c, m) in centres.iter().zip(&members) {
            match (m.len() > capacity).then(|| split(m)).flatten() {
                Some(pair) => {
                    next.extend(pair);
                    split_any = true;
                }
                None => next.push(*c),
            }
        }
        if !split_any {
            let load = members.iter().map(Vec::len).collect();
            return Ok((cs, Report { rounds, load }));
        }
        centres = next;
    }
    unreachable!()
}
