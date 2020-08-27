//! A test model with energy linearly dependent on a single parameter.

use super::*;

use rand::prelude::*;

/// The parameters needed to configure an Linear model.
///
/// These parameters are normally set via command-line arguments.
#[derive(Serialize, Deserialize, Debug, AutoArgs)]
#[allow(non_snake_case)]
pub struct LinearParams {
    /// Max distance
    pub L: Length,
}

#[allow(non_snake_case)]
/// An Linear model.
#[derive(Serialize, Deserialize, Debug)]
pub struct Linear {
    /// The energy of the system
    pub x: Length,
    /// The dimensions of the box.
    pub L: Length,
    /// The last change we made (and might want to undo).
    possible_change: Option<Length>,
}

impl From<LinearParams> for Linear {
    fn from(params: LinearParams) -> Linear {
        Linear {
            x: 0.5*params.L,
            L: params.L,
            possible_change: None,
        }
    }
}

impl System for Linear {
    fn energy(&self) -> Energy {
        self.x/Length::new(1.)*Energy::new(1.)
    }
    fn compute_energy(&self) -> Energy {
        self.energy()
    }
}

impl ConfirmSystem for Linear {
    fn confirm(&mut self) {
        if let Some(x) = self.possible_change {
            self.x = x;
        }
    }
}

impl MovableSystem for Linear {
    fn plan_move(&mut self, rng: &mut MyRng, d: Length) -> Option<Energy> {
        let v: f64 = rng.sample(rand_distr::StandardNormal);
        let xnew = self.x + d*v;
        if xnew > Length::new(0.) && xnew < self.L {
            self.possible_change = Some(xnew);
            Some(xnew*Energy::new(1.)/Length::new(1.))
        } else {
            None
        }
    }
}

