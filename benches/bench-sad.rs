#[macro_use]
extern crate criterion;
extern crate sadmc;
extern crate rand;

use rand::Rng;
use criterion::Criterion;

fn criterion_benchmark(c: &mut Criterion) {
    let mut rng = sadmc::rng::MyRng::from_u64(0);
    c.bench_function("MyRng.gen<u64>", move |b| b.iter(|| rng.gen::<u64>()));
    let mut rng = sadmc::rng::MyRng::from_u64(0);
    c.bench_function("MyRng.gen<f64>", move |b| b.iter(|| rng.gen::<f64>()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
