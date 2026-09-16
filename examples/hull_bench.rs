//! Time both Delaunay backends on random and Fibonacci-lattice sites.
//! Usage: cargo run --release --example hull_bench

use hipparchus_rs::{V3, build::fibonacci, hull};
use std::time::Instant;

fn main() {
    let mut s = 0x2545F4914F6CDD1Du64;
    let mut f = move || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    for n in [1_000, 5_000, 20_000, 100_000, 1_000_000] {
        let pts: Vec<V3> = (0..n)
            .map(|_| {
                let z = 2.0 * f() - 1.0;
                let t = std::f64::consts::TAU * f();
                let r = (1.0 - z * z).sqrt();
                V3::new(r * t.cos(), r * t.sin(), z)
            })
            .collect();
        let t = Instant::now();
        hull::sphere_delaunay_geo(&pts).unwrap();
        let g = t.elapsed();
        let fib = fibonacci(n);
        let t = Instant::now();
        hull::sphere_delaunay_geo(&fib).unwrap();
        let gf = t.elapsed();
        let r = if n <= 20_000 {
            let t = Instant::now();
            hull::sphere_delaunay(&pts).unwrap();
            format!("{:.2?}", t.elapsed())
        } else {
            "-".into()
        };
        println!("n={n:>9}  robust {r:>10}  geo {g:>10.2?}  geo(fibonacci) {gf:>10.2?}");
    }
}
