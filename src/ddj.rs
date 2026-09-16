//! Readers for the DDJ 1992 `COAST*.IDX` / `COAST.DAT` files.

use crate::algebra::{Frame, V3};

/// Local units per object radius: `value = plane × LC_SCALE / radius`.
pub const LC_SCALE: f64 = 32000.0;

/// 40-byte little-endian index record (16-bit MS C 6.0 layout).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IndexRec {
    /// Object centre.
    pub centre: V3,
    /// Object radius in radians.
    pub radius: f64,
    /// Byte offset of the first vertex in `COAST.DAT`.
    pub offset: u32,
    /// Number of vertices.
    pub count: u16,
    /// Segment identifier.
    pub id: i16,
}

fn f64_at(b: &[u8], i: usize) -> f64 {
    f64::from_le_bytes(b[i..i + 8].try_into().unwrap())
}

/// Reads the records of a `COAST*.IDX` file. Bytes after the last complete
/// record are ignored.
pub fn read_index(bytes: &[u8]) -> Vec<IndexRec> {
    bytes
        .as_chunks::<40>()
        .0
        .iter()
        .map(|r| IndexRec {
            centre: V3::new(f64_at(r, 0), f64_at(r, 8), f64_at(r, 16)),
            radius: f64_at(r, 24),
            offset: u32::from_le_bytes(r[32..36].try_into().unwrap()),
            count: u16::from_le_bytes(r[36..38].try_into().unwrap()),
            id: i16::from_le_bytes(r[38..40].try_into().unwrap()),
        })
        .collect()
}

/// Decode one segment's `(i16, i16)` local vertices, as `VIEW.C` does.
///
/// # Panics
///
/// If `dat` does not hold `rec.count` vertices from `rec.offset`.
pub fn decode_segment(rec: &IndexRec, dat: &[u8]) -> Vec<V3> {
    let frame = Frame::new(rec.centre);
    let s = rec.radius / LC_SCALE;
    let start = rec.offset as usize;
    dat[start..start + 4 * rec.count as usize]
        .as_chunks::<4>()
        .0
        .iter()
        .map(|v| {
            let e = i16::from_le_bytes([v[0], v[1]]) as f64;
            let n = i16::from_le_bytes([v[2], v[3]]) as f64;
            frame.unmap([s * e, s * n])
        })
        .collect()
}
