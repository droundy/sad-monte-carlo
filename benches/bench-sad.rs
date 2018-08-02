#[macro_use]
extern crate criterion;
extern crate sadmc;
extern crate rand;

use rand::Rng;
use criterion::Criterion;

use sadmc::mc::MonteCarlo;
use sadmc::mc::sad::{Sad, SadParams};
use sadmc::system::square::{SquareWell, SquareWellNParams};
use sadmc::system::{Length, MovableSystem};

fn gen_sw(n_atoms: usize) -> SquareWell {
    let mut sw_params = SquareWellNParams::default();
    sw_params.N = n_atoms;
    let mut sw = SquareWell::from(sw_params);
    // Randomize things a bit before beginning.
    let mut rng = sadmc::rng::MyRng::from_u64(1);
    for _ in 0..n_atoms*1000 {
        sw.move_once(&mut rng, Length::new(1.0));
    }
    sw
}

fn gen_sad(n_atoms: usize) -> Sad<SquareWell> {
    let sad_params = SadParams::default();
    let mut sad = Sad::<SquareWell>::from_params(sad_params, gen_sw(n_atoms),
                                                 ::std::path::PathBuf::from("/dev/null"));
    // Randomize things a bit before beginning.
    for _ in 0..n_atoms*1000 {
        sad.move_once();
    }
    sad
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut rng = sadmc::rng::MyRng::from_u64(0);
    c.bench_function("MyRng.gen<u64>", move |b| b.iter(|| rng.gen::<u64>()));
    let mut rng = sadmc::rng::MyRng::from_u64(0);
    c.bench_function("MyRng.gen<f64>", move |b| b.iter(|| rng.gen::<f64>()));

    c.bench_function_over_inputs("sad_sw_move_once",
                                 move |b, &&n_atoms| {
                                     let mut sad = gen_sad(n_atoms);
                                     b.iter(|| sad.move_once())
                                 },
                                 &[50, 100,
                                   // 200, 400
                                   ]);

    c.bench_function_over_inputs("SW_sw_move_once",
                                 move |b, &&n_atoms| {
                                     let mut sw = gen_sw(n_atoms);
                                     let mut rng = sadmc::rng::MyRng::from_u64(2);
                                     b.iter(|| sw.move_once(&mut rng, Length::new(0.1)))
                                 },
                                 &[50, 100,
                                   // 200, 400
                                   ]);
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
