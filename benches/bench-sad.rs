#[macro_use]
extern crate criterion;
extern crate sadmc;
extern crate rand;

use rand::Rng;
use rand::distributions::Uniform;
use criterion::Criterion;

use sadmc::system::units;
use sadmc::mc::MonteCarlo;
use sadmc::mc::sad::{Sad, SadParams};
use sadmc::mc::energy::{EnergyMC, EnergyMCParams};
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

fn gen_energy_sad(n_atoms: usize) -> EnergyMC<SquareWell> {
    let params = EnergyMCParams::default();
    let mut mc = EnergyMC::<SquareWell>::from_params(params, gen_sw(n_atoms),
                                                     ::std::path::PathBuf::from("/dev/null"));
    // Randomize things a bit before beginning.
    for _ in 0..n_atoms*1000 {
        mc.move_once();
    }
    mc
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

    c.bench_function_over_inputs("energy_sad_sw_move_once",
                                 move |b, &&n_atoms| {
                                     let mut sad = gen_energy_sad(n_atoms);
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

    let closest = criterion::Fun::new("standard", |b,&n_atoms| {
        let sw = gen_sw(n_atoms);
        let mut rng = sadmc::rng::MyRng::from_u64(2);
        b.iter_with_setup(|| {
            let r1 = sw.positions[rng.sample(Uniform::new(0, sw.positions.len()))];
            let r2 = sw.positions[rng.sample(Uniform::new(0, sw.positions.len()))];
            (r1,r2)
        }, |(r1,r2)| sw.closest_distance2(r1,r2));
    });
    let unsafe_closest = criterion::Fun::new("unsafe", |b,&n_atoms| {
        let sw = gen_sw(n_atoms);
        let mut rng = sadmc::rng::MyRng::from_u64(2);
        b.iter_with_setup(|| {
            let r1 = sw.positions[rng.sample(Uniform::new(0, sw.positions.len()))];
            let r2 = sw.positions[rng.sample(Uniform::new(0, sw.positions.len()))];
            (r1,r2)
        }, |(r1,r2)| sw.unsafe_closest_distance2(r1,r2));
    });
    let sloppy_closest = criterion::Fun::new("sloppy", |b,&n_atoms| {
        let sw = gen_sw(n_atoms);
        let mut rng = sadmc::rng::MyRng::from_u64(2);
        b.iter_with_setup(|| {
            let r1 = sw.positions[rng.sample(Uniform::new(0, sw.positions.len()))];
            let r2 = sw.positions[rng.sample(Uniform::new(0, sw.positions.len()))];
            (r1,r2)
        }, |(r1,r2)| sw.sloppy_closest_distance2(r1,r2));
    });
    let funs = vec![closest, unsafe_closest, sloppy_closest];
    c.bench_functions("SW_closest_distance2", funs, 50);

    let put = criterion::Fun::new("standard", |b,&n_atoms| {
        let sw = gen_sw(n_atoms);
        let mut rng = sadmc::rng::MyRng::from_u64(2);
        b.iter_with_setup(|| {
            let r = sw.positions[rng.sample(Uniform::new(0, sw.positions.len()))];
            r + rng.vector()*0.1*units::SIGMA
        }, |r| sw.put_in_cell(r));
    });
    let sloppy_put = criterion::Fun::new("sloppy", |b,&n_atoms| {
        let sw = gen_sw(n_atoms);
        let mut rng = sadmc::rng::MyRng::from_u64(2);
        b.iter_with_setup(|| {
            let r = sw.positions[rng.sample(Uniform::new(0, sw.positions.len()))];
            r + rng.vector()*0.1*units::SIGMA
        }, |r| sw.sloppy_put_in_cell(r));
    });
    let funs = vec![put, sloppy_put];
    c.bench_functions("SW_put_in_cell", funs, 50);
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
