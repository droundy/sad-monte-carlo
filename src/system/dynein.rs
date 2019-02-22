//! A dynein walk.

use super::*;

use rand::prelude::*;

/// The parameters needed to configure dynein walking.
#[derive(Serialize, Deserialize, Debug, ClapMe)]
#[allow(non_snake_case)]
pub struct DyneinParams {
    /// Motor Domain Angles
    nm: usize,
    fm: usize,
    ///Initial Displacement
    initd: usize,

}

#[allow(non_snake_case)]
/// A dynein model.
#[derive(Serialize, Deserialize, Debug)]
pub struct Dynein {
    /// The energy of the system
    E: Energy,
    /// The unbinding rates
    K: f64,
    /// The motor domain angles
    pub nm: usize,
    pub fm: usize,
    /// The initial displacement of foot domain
    pub initd: Length,
    /// The stalk length.
    Ls: Length,  
    /// The tail length.
    Lt: Length, 
    /// The spins themselves (Change This)
    ///S: Vec<i8>, (Disregard this)   
    /// The unbinding constant
    kub: usize,
    /// The equilbrium Angles
    eqb: usize,
    eqmpre: usize,
    eqmpost: usize,
    eqt: usize,


    /// The last change we made (and might want to undo).
    possible_change: Option<(usize, Energy)>,
}

impl From<DyneinParams> for Dynein {
    fn from(params: DyneinParams) -> Dynein {
        let mut dynein = Dynein {
            E: Energy::new(0.),
            nm: params.nm,
            fm: params.fm,
            initd: params.initd,
           //S: vec![1; params.N*params.N],
            Ls = 20.75,     //nm
            Lt = 23,        //nm
            kub = 30,       //s^-1
            eqb = 120,      
            eqmpre = 197,
            eqmpost = 242,
            eqt = 0,

            possible_change: None,
        };
        //assert!(dynein.N > 1); // otherwise, we are our own neighbor!
        let mut rng = ::rng::MyRng::from_u64(10137);
        /*for s in dynein.S.iter_mut() {
            *s = *rng.choose(&[-1,1]).unwrap();
        }*/
        dynein.E = dynein.compute_energy();
        dynein
    }
}

impl Dynien {
    fn unbinding_rate(%self) -> f64 {
        self.kub * E.powf(self.b*(self.nm-self.eqmpre)) // This is wrong, b = 1/(Boltzmann Constant * T) ?
    }

    fn  
}

impl System for Dynein {
    fn energy(&self) -> Energy {
        self.E
    }
    fn compute_energy(&self) -> Energy {
        let mut e: Energy = units::EPSILON*0.0;
        /*
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
        }*/
        e
    }
    fn delta_energy(&self) -> Option<Energy> {
        Some(4.*units::EPSILON)
    }
}

impl ConfirmSystem for Dynein {
    fn confirm(&mut self) {
        if let Some((i, e)) = self.possible_change {
            self.S[i] *= -1;
            self.E = e;
        }
    }
}

impl MovableSystem for Dynein {
    fn plan_move(&mut self, rng: &mut MyRng, _: Length) -> Option<Energy> {
        /*
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
        */
    }
}

#[cfg(test)]
#[allow(non_snake_case)]
fn energy_works_with_N(N: usize) {
    let mut dynein = dynein::from(DyneinParams { N: N });

    println!("starting energy...");
    assert_eq!(dynein.energy(), dynein.compute_energy());

    let mut rng = ::rng::MyRng::from_u64(10137);
    for _ in 0..10000 {
        dynein.plan_move(&mut rng, Length::new(0.0));
        dynein.confirm();
        assert_eq!(dynein.energy(), dynein.compute_energy());
    }
}

#[test]
fn energy_works() {
    for &n in &[2,3,10,15,137,150] {
        println!("testing with N={}", n);
        energy_works_with_N(n);
    }
}

