//! A square well fluid.

use super::*;

use rand::prelude::*;

/// The parameters needed to configure a square well system.
#[derive(Serialize, Deserialize, Debug, ClapMe)]
#[allow(non_snake_case)]
pub struct IsingParams {
    /// Width of the square grid
    N: usize,
}

#[allow(non_snake_case)]
/// A square well fluid.
#[derive(Serialize, Deserialize, Debug)]
pub struct Ising {
    /// The energy of the system
    E: Energy,
    /// The dimensions of the box.
    pub N: usize,
    /// The spins themselves
    S: Vec<i8>,
    /// The last change we made (and might want to undo).
    last_changed: Option<(usize, Energy)>,
}

impl From<IsingParams> for Ising {
    fn from(params: IsingParams) -> Ising {
        let mut ising = Ising {
            E: Energy::new(0.),
            N: params.N,
            S: vec![1; params.N*params.N],
            last_changed: None,
        };
        assert!(ising.N > 1); // otherwise, we are our own neighbor!
        let mut rng = ::rng::MyRng::from_u64(10137);
        for s in ising.S.iter_mut() {
            *s = *rng.choose(&[-1,1]).unwrap();
        }
        ising.E = ising.compute_energy();
        ising
    }
}

impl System for Ising {
    fn energy(&self) -> Energy {
        self.E
    }
    fn compute_energy(&self) -> Energy {
        let mut e: Energy = units::EPSILON*0.0;
        for i1 in 0 .. self.N {
            for j1 in 0 .. self.N {
                let j2 = (j1 + 1) % self.N;
                let mut neighbor_tot = self.S[i1 + j2*self.N];
                // let j2 = modulus(j - 1, self.N);
                // neighbor_tot += self.S[i1 + j2*self.N];
                let i2 = (i1 + 1) % self.N;
                neighbor_tot += self.S[i2 + j1*self.N];
                // let i2 = modulus(i - 1, self.N);
                // neighbor_tot += self.S[i1 + j2*self.N];
                e += neighbor_tot as f64 * Energy::new(self.S[i1+j1*self.N] as f64);
            }
        }
        e
    }
    fn delta_energy(&self) -> Option<Energy> {
        Some(4.*units::EPSILON)
    }
}

impl UndoSystem for Ising {
    fn undo(&mut self) {
        if let Some((i, e)) = self.last_changed {
            self.S[i] *= -1;
            self.E = e;
        }
    }
}

impl MovableSystem for Ising {
    fn move_once(&mut self, rng: &mut MyRng, _: Length) -> Option<Energy> {
        let i = rng.gen_range(0, self.N);
        let j = rng.gen_range(0, self.N);
        self.S[i+j*self.N] *= -1;
        self.last_changed = Some((i+j*self.N, self.E));

        let j2 = (j + 1) % self.N;
        let mut neighbor_tot = self.S[i + j2*self.N];
        let j2 = (j + self.N - 1) % self.N;
        neighbor_tot += self.S[i + j2*self.N];
        let i2 = (i + 1) % self.N;
        neighbor_tot += self.S[i2 + j*self.N];
        let i2 = (i + self.N - 1) % self.N;
        neighbor_tot += self.S[i2 + j*self.N];
        let e = neighbor_tot as f64 * self.S[i+j*self.N] as f64 * Energy::new(2.0);
        self.E += e;
        Some(e)
    }
}

#[cfg(test)]
#[allow(non_snake_case)]
fn energy_works_with_N(N: usize) {
    let mut ising = Ising::from(IsingParams { N: N });

    println!("starting energy...");
    assert_eq!(ising.energy(), ising.compute_energy());

    let mut rng = ::rng::MyRng::from_u64(10137);
    for _ in 0..10000 {
        ising.move_once(&mut rng, Length::new(0.0));
        assert_eq!(ising.energy(), ising.compute_energy());
    }
}

#[test]
fn energy_works() {
    for &n in &[2,3,10,15,137,150] {
        println!("testing with N={}", n);
        energy_works_with_N(n);
    }
}
