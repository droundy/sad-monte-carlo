//! A dynein walk.

use super::*;

use rand::prelude::*;

/// The parameters needed to configure a dynein system.
#[derive(Serialize, Deserialize, Debug, ClapMe)]
#[allow(non_snake_case)]
pub struct DyneinParams {
    /// Stalk length in nm
    pub Ls: f64,
    /// Tail length in nm
    pub Lt: f64,
}

impl Default for DyneinParams {
    fn default() -> Self {
        DyneinParams {
            Ls: 11.0,
            Lt: 23.0,
        }
    }
}

#[allow(non_snake_case)]
/// A square well fluid.
#[derive(Serialize, Deserialize, Debug)]
pub struct Dynein {
    /// The parameters defining our model
    pub params: DyneinParams,
    /// The energy of the system
    E: Energy,
    /// The last change we made (and might want to undo).
    possible_change: Option<(usize, Energy)>,
}

impl From<DyneinParams> for Dynein {
    fn from(params: DyneinParams) -> Dynein {
        let mut dynein = Dynein {
            params,
            E: Energy::new(0.),
            possible_change: None,
        };
        dynein.E = dynein.compute_energy();
        dynein
    }
}

impl System for Dynein {
    fn energy(&self) -> Energy {
        self.E
    }
    fn compute_energy(&self) -> Energy {
        let mut e: Energy = units::EPSILON*0.0;
        // FIXME need to compute the energy here
        e
    }
}

impl ConfirmSystem for Dynein {
    fn confirm(&mut self) {
        if let Some((i, e)) = self.possible_change {
            // FIXME need to make a change...
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
        None // FIXME need to compute actual move and energy change
    }
}

#[test]
fn energy_works() {
    let mut dynein = Dynein::from(DyneinParams::default());

    println!("starting energy...");
    assert_eq!(dynein.energy(), dynein.compute_energy());

    // let mut rng = ::rng::MyRng::from_u64(10137);
    // for _ in 0..10000 {
    //     Dynein.plan_move(&mut rng, Length::new(0.0));
    //     Dynein.confirm();
    //     assert_eq!(Dynein.energy(), Dynein.compute_energy());
    // }
}
