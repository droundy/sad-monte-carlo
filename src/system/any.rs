//! A "system" that can represent any system

use super::*;

/// The parameters needed to configure an Any model.
///
/// These parameters are normally set via command-line arguments.
#[derive(Serialize, Deserialize, Debug, AutoArgs)]
#[allow(non_snake_case)]
pub enum AnyParams {
    /// a fake system
    Fake(fake::FakeParams),
    /// a fake erfinv system
    FakeErfinv(erfinv::ParametersN),
    /// a wca system
    Wca(wca::WcaNParams),
    /// a lj system
    Lj(lj::LjParams),
    /// an Ising system
    Ising(ising::IsingParams),
    /// a square well system
    Sw(optsquare::SquareWellNParams),
}

#[allow(non_snake_case)]
/// An Any model.
#[derive(Serialize, Deserialize, Debug)]
pub enum Any {
    /// A fake system
    Fake(fake::Fake),
    /// A fake erfinv system
    FakeErfinv(erfinv::ErfInv),
    /// A wca system
    Wca(wca::Wca),
    /// A lj system
    Lj(lj::Lj),
    /// An Ising system
    Ising(ising::Ising),
    /// a square well system
    Sw(optsquare::SquareWell),
}

impl From<AnyParams> for Any {
    fn from(parameters: AnyParams) -> Any {
        match parameters {
            AnyParams::Fake(parameters) => Any::Fake(fake::Fake::from(parameters)),
            AnyParams::Wca(parameters) => Any::Wca(wca::Wca::from(parameters)),
            AnyParams::Lj(parameters) => Any::Lj(lj::Lj::from(parameters)),
            AnyParams::Ising(parameters) => Any::Ising(ising::Ising::from(parameters)),
            AnyParams::FakeErfinv(parameters) => Any::FakeErfinv(erfinv::ErfInv::from(parameters)),
            AnyParams::Sw(parameters) => Any::Sw(optsquare::SquareWell::from(parameters)),
        }
    }
}

impl Any {
    fn movable(&self) -> &dyn MovableSystem {
        match self {
            Any::Fake(s) => s as &dyn MovableSystem,
            Any::Wca(s) => s as &dyn MovableSystem,
            Any::Lj(s) => s as &dyn MovableSystem,
            Any::Ising(s) => s as &dyn MovableSystem,
            Any::FakeErfinv(s) => s as &dyn MovableSystem,
            Any::Sw(s) => s as &dyn MovableSystem,
        }
    }
    fn movable_mut(&mut self) -> &mut dyn MovableSystem {
        match self {
            Any::Fake(s) => s as &mut dyn MovableSystem,
            Any::Wca(s) => s as &mut dyn MovableSystem,
            Any::Lj(s) => s as &mut dyn MovableSystem,
            Any::Ising(s) => s as &mut dyn MovableSystem,
            Any::FakeErfinv(s) => s as &mut dyn MovableSystem,
            Any::Sw(s) => s as &mut dyn MovableSystem,
        }
    }
}

impl System for Any {
    fn energy(&self) -> Energy {
        self.movable().energy()
    }
    fn compute_energy(&self) -> Energy {
        self.movable().compute_energy()
    }
}

impl ConfirmSystem for Any {
    fn confirm(&mut self) {
        self.movable_mut().confirm()
    }
}

impl MovableSystem for Any {
    fn plan_move(&mut self, rng: &mut MyRng, d: Length) -> Option<Energy> {
        self.movable_mut().plan_move(rng, d)
    }
}
