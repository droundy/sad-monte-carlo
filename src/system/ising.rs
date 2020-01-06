//! The Ising model

use super::*;

use rand::prelude::*;

/// The parameters needed to configure an Ising model.
///
/// These parameters are normally set via command-line arguments.
#[derive(Serialize, Deserialize, Debug, ClapMe)]
#[allow(non_snake_case)]
pub struct IsingParams {
    /// Width of the square grid
    pub N: usize,
}

#[allow(non_snake_case)]
/// An Ising model.
#[derive(Serialize, Deserialize, Debug)]
pub struct Ising {
    /// The energy of the system
    E: Energy,
    /// The dimensions of the box.
    pub N: usize,
    /// The spins themselves
    S: Vec<i8>,
    /// The last change we made (and might want to undo).
    possible_change: Option<(usize, Energy)>,
}

impl From<IsingParams> for Ising {
    fn from(params: IsingParams) -> Ising {
        let mut ising = Ising {
            E: Energy::new(0.),
            N: params.N,
            S: vec![1; params.N*params.N],
            possible_change: None,
        };
        assert!(ising.N > 1); // otherwise, we are our own neighbor!
        // This is a bit hokey.  We have a fixed random number seed
        // for initializing the spins randomly.  This is because this
        // function doesn't have access to the seed used by the MC
        // simulation, and it seems excessive to have two distinct
        // seeds.
        let mut rng = crate::rng::MyRng::from_u64(10137);
        for s in ising.S.iter_mut() {
            // The following picks a random number +1 or -1
            *s = (rng.next_u64() as i8 & 1)*2 - 1;
        }
        ising.E = ising.compute_energy();
        ising
    }
}

impl System for Ising {
    type CollectedData = ();
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

impl ConfirmSystem for Ising {
    fn confirm(&mut self) {
        if let Some((i, e)) = self.possible_change {
            self.S[i] *= -1;
            self.E = e;
        }
    }
}

impl MovableSystem for Ising {
    fn plan_move(&mut self, rng: &mut MyRng, _: Length) -> Option<Energy> {
        let i = rng.gen_range(0, self.N);
        let j = rng.gen_range(0, self.N);

        let j2 = (j + 1) % self.N;
        let mut neighbor_tot = self.S[i + j2*self.N];
        let j2 = (j + self.N - 1) % self.N;
        neighbor_tot += self.S[i + j2*self.N];
        let i2 = (i + 1) % self.N;
        neighbor_tot += self.S[i2 + j*self.N];
        let i2 = (i + self.N - 1) % self.N;
        neighbor_tot += self.S[i2 + j*self.N];
        let e = self.E - neighbor_tot as f64 * self.S[i+j*self.N] as f64 * Energy::new(2.0);
        self.possible_change = Some((i+j*self.N, e));
        Some(e)
    }
}

#[cfg(test)]
#[allow(non_snake_case)]
fn energy_works_with_N(N: usize) {
    let mut ising = Ising::from(IsingParams { N: N });

    println!("starting energy...");
    assert_eq!(ising.energy(), ising.compute_energy());

    let mut rng = crate::rng::MyRng::from_u64(10137);
    for _ in 0..10000 {
        ising.plan_move(&mut rng, Length::new(0.0));
        ising.confirm();
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
