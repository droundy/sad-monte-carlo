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
    /// a water system
    Water(water::WaterParams),
    /// an Ising system
    Ising(ising::IsingParams),
    /// a square well system
    Sw(optsquare::SquareWellNParams),
}
/// The parameters needed to configure an AnyGrand model.
///
/// These parameters are normally set via command-line arguments.
#[derive(Serialize, Deserialize, Debug, AutoArgs)]
#[allow(non_snake_case)]
pub enum AnyGrandParams {
    /// a fake erfinv system
    FakeErfinv(erfinv::Parameters),
    /// a wca system
    Wca(wca::WcaParams),
    /// a lj system
    Lj(lj::GrandLjParams),
    /// a square well system
    Sw(optsquare::SquareWellParams),
}

#[allow(non_snake_case)]
/// An AnyGrand model which can be used with Grand MC.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum AnyGrand {
    /// A fake erfinv system
    FakeErfinv(erfinv::ErfInv),
    /// A wca system
    Wca(wca::Wca),
    /// A lj system
    Lj(lj::Lj),
    /// a square well system
    Sw(optsquare::SquareWell),
}

#[allow(non_snake_case)]
/// An Any model.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum Any {
    /// A fake system
    Fake(fake::Fake),
    /// A fake erfinv system
    FakeErfinv(erfinv::ErfInv),
    /// A wca system
    Wca(wca::Wca),
    /// A lj system
    Lj(lj::Lj),
    /// A water system
    Water(water::Water),
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
            AnyParams::Water(parameters) => Any::Water(water::Water::from(parameters)),
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
            Any::Water(s) => s as &dyn MovableSystem,
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
            Any::Water(s) => s as &mut dyn MovableSystem,
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
    fn randomize(&mut self, rng: &mut MyRng) -> Energy {
        self.movable_mut().randomize(rng)
    }
    fn update_caches(&mut self) {
        self.movable_mut().update_caches();
    }
    fn min_moves_to_randomize(&self) -> u64 {
        self.movable().min_moves_to_randomize()
    }
    fn dimensionality(&self) -> u64 {
        self.movable().dimensionality()
    }
    fn data_to_collect(&self, iter: u64) -> Vec<(Interned, f64)> {
        self.movable().data_to_collect(iter)
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
    fn max_size(&self) -> Length {
        self.movable().max_size()
    }
}

impl From<AnyGrandParams> for AnyGrand {
    fn from(parameters: AnyGrandParams) -> AnyGrand {
        match parameters {
            AnyGrandParams::Wca(parameters) => AnyGrand::Wca(wca::Wca::from(parameters)),
            AnyGrandParams::Lj(parameters) => AnyGrand::Lj(lj::Lj::from(parameters)),
            AnyGrandParams::FakeErfinv(parameters) => {
                AnyGrand::FakeErfinv(erfinv::ErfInv::from(parameters))
            }
            AnyGrandParams::Sw(parameters) => AnyGrand::Sw(optsquare::SquareWell::from(parameters)),
        }
    }
}

impl AnyGrand {
    fn grand(&self) -> &dyn GrandSystem {
        match self {
            AnyGrand::Wca(s) => s as &dyn GrandSystem,
            AnyGrand::Lj(s) => s as &dyn GrandSystem,
            AnyGrand::FakeErfinv(s) => s as &dyn GrandSystem,
            AnyGrand::Sw(s) => s as &dyn GrandSystem,
        }
    }
    fn grand_mut(&mut self) -> &mut dyn GrandSystem {
        match self {
            AnyGrand::Wca(s) => s as &mut dyn GrandSystem,
            AnyGrand::Lj(s) => s as &mut dyn GrandSystem,
            AnyGrand::FakeErfinv(s) => s as &mut dyn GrandSystem,
            AnyGrand::Sw(s) => s as &mut dyn GrandSystem,
        }
    }
}

impl System for AnyGrand {
    fn energy(&self) -> Energy {
        self.grand().energy()
    }
    fn compute_energy(&self) -> Energy {
        self.grand().compute_energy()
    }
    fn randomize(&mut self, rng: &mut MyRng) -> Energy {
        self.grand_mut().randomize(rng)
    }
    fn min_moves_to_randomize(&self) -> u64 {
        self.grand().min_moves_to_randomize()
    }
    fn dimensionality(&self) -> u64 {
        self.grand().dimensionality()
    }
}

impl ConfirmSystem for AnyGrand {
    fn confirm(&mut self) {
        self.grand_mut().confirm()
    }
}

impl MovableSystem for AnyGrand {
    fn plan_move(&mut self, rng: &mut MyRng, d: Length) -> Option<Energy> {
        self.grand_mut().plan_move(rng, d)
    }
    fn max_size(&self) -> Length {
        self.grand().max_size()
    }
}

impl GrandSystem for AnyGrand {
    fn plan_add(&mut self, rng: &mut MyRng) -> Option<Energy> {
        self.grand_mut().plan_add(rng)
    }

    fn plan_remove(&mut self, rng: &mut MyRng) -> Energy {
        self.grand_mut().plan_remove(rng)
    }

    fn num_atoms(&self) -> usize {
        self.grand().num_atoms()
    }
}

impl GrandReplicaSystem for AnyGrand {
    fn plan_swap_atom(&self, other: &Self, rng: &mut MyRng) -> Option<(usize, Energy, Energy)> {
        use AnyGrand::*;
        match (self, other) {
            (Wca(s), Wca(o)) => s.plan_swap_atom(o, rng),
            (Lj(s), Lj(o)) => s.plan_swap_atom(o, rng),
            (FakeErfinv(s), FakeErfinv(o)) => s.plan_swap_atom(o, rng),
            (Sw(s), Sw(o)) => s.plan_swap_atom(o, rng),
            _ => unreachable!(),
        }
    }

    fn swap_atom(&mut self, other: &mut Self, which: usize) {
        use AnyGrand::*;
        match (self, other) {
            (Wca(s), Wca(o)) => s.swap_atom(o, which),
            (Lj(s), Lj(o)) => s.swap_atom(o, which),
            (FakeErfinv(s), FakeErfinv(o)) => s.swap_atom(o, which),
            (Sw(s), Sw(o)) => s.swap_atom(o, which),
            _ => unreachable!(),
        }
    }
}
