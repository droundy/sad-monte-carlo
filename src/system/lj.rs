//! A lennard jones fluid.

use super::*;

use crate::prettyfloat::PrettyFloat;
use dimensioned::{Abs, Dimensionless};
use rand::distributions::Uniform;
use rand::{Rng, SeedableRng};
use vector3d::Vector3d;

/// The parameters needed to configure a lennard jones  system.
#[allow(non_snake_case)]
#[derive(Serialize, Deserialize, Debug, AutoArgs)]
pub struct LjParams {
    /// The number of atoms
    N: usize,
    /// The radius of the spherical box
    radius: Length,
}

/// The parameters needed to configure a lennard jones system for a grand computation.
#[allow(non_snake_case)]
#[derive(Serialize, Deserialize, Debug, AutoArgs)]
pub struct GrandLjParams {
    /// The radius of the spherical box
    radius: Length,
}

#[allow(non_snake_case)]
/// A lennard jones fluid.
#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct Lj {
    /// The energy of the system
    E: Energy,
    /// The estimated accumulated error so far in E
    error: Energy,
    /// The last change we made (and might want to undo).
    possible_change: Change,
    /// The locations of the atoms
    positions: Vec<Vector3d<Length>>,
    /// The square of the maximum radius permitted
    max_radius_squared: Area,
    /// The maximum radius permitted
    max_radius: Length,
}

/// This defines the energy/radial bins.
#[derive(Serialize, Deserialize, Debug, Default,Clone)]
pub struct Collected {
    /// The number of atoms for each radial bin from center of sphere.
    pub from_center: Vec<u64>,
    /// The number of atoms for each radial bin from center of masse.
    pub from_cm: Vec<u64>,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
/// Define the types of changes that can be made to the system
enum Change {
    /// Move an atom already in the system
    Move {
        which: usize,
        to: Vector3d<Length>,
        e: Energy,
    },
    Add {
        to: Vector3d<Length>,
        e: Energy,
    },
    Remove {
        which: usize,
        e: Energy,
    },
    /// Make no changes to the system
    None,
}

/// Define the WCA interaction potential and criteria
fn potential(r_squared: Area) -> Energy {
    let sig_sqr = units::SIGMA * units::SIGMA;
    4.0 * units::EPSILON * ((sig_sqr / r_squared).powi(6) - (sig_sqr / r_squared).powi(3))
}

impl Lj {
    /// Move a specified atom.  Returns the change in energy, or
    /// `None` if the atom could not be placed there.
    pub fn move_atom(&mut self, which: usize, r: Vector3d<Length>) -> Option<Energy> {
        let previous_rsqr = self.positions[which].norm2();
        if r.norm2() > self.max_radius_squared && r.norm2() > previous_rsqr {
            return None;
        }
        let mut e = self.E;
        let from = self.positions[which];
        for r1 in self
            .positions
            .iter()
            .cloned()
            .enumerate()
            .filter(|&(i, _)| i != which)
            .map(|(_, x)| x)
        {
            e += potential((r1 - r).norm2()) - potential((r1 - from).norm2());
        }
        self.possible_change = Change::Move { which, to: r, e };
        Some(e)
    }
    fn expected_accuracy(&self, newe: Energy) -> Energy {
        newe.abs() * 1e-14 * (self.positions.len() as f64) * (self.positions.len() as f64)
    }

    fn set_energy(&mut self, new_e: Energy) {
        let new_error = if new_e.abs() > self.E.abs() {
            new_e.abs() * 1e-15 * (self.positions.len() as f64)
        } else {
            self.E.abs() * 1e-15 * (self.positions.len() as f64)
        };
        self.error = new_error + self.error;
        if self.error > self.expected_accuracy(new_e) {
            self.error *= 0.0;
            self.E = self.compute_energy();
        } else {
            self.E = new_e;
        }
    }
}

impl From<LjParams> for Lj {
    fn from(params: LjParams) -> Lj {
        let mut rng = crate::rng::MyRng::seed_from_u64(0);
        let mut best_energy = 1e80 * units::EPSILON;
        let mut best_positions = Vec::new();
        println!("I am creating a LJ system with {} atoms!", params.N);
        for attempt in 0..10000000 {
            let mut positions = Vec::new();
            for _ in 0..params.N {
                let mut r;
                loop {
                    r = Vector3d::new(
                        rng.sample(Uniform::new(-1.0, 1.0)),
                        rng.sample(Uniform::new(-1.0, 1.0)),
                        rng.sample(Uniform::new(-1.0, 1.0)),
                    );
                    if r.norm2() < 1.0 {
                        break;
                    }
                }
                positions.push(r * params.radius);
            }
            let mut cm = positions[0];
            for x in positions.iter().skip(1) {
                cm = cm + *x;
            }
            cm = cm / params.N as f64;
            for x in positions.iter_mut() {
                *x = *x - cm;
            }
            if positions
                .iter()
                .any(|x| x.norm2() > params.radius * params.radius)
            {
                continue;
            }
            let mut lj = Lj {
                E: 0.0 * units::EPSILON,
                error: 0.0 * units::EPSILON,
                possible_change: Change::None,
                positions,
                max_radius_squared: params.radius * params.radius,
                max_radius: params.radius,
            };
            lj.E = lj.compute_energy();
            if lj.E < best_energy {
                best_energy = lj.E;
                best_positions = lj.positions.clone();
                println!(
                    "found a new best energy{:?} after {} attempts",
                    lj.E, attempt
                );
            }
            if lj.E < 0.0 * units::EPSILON {
                return lj;
            }
        }
        let mut lj = Lj {
            E: best_energy,
            error: 0.0 * units::EPSILON,
            possible_change: Change::None,
            positions: best_positions,
            max_radius_squared: params.radius * params.radius,
            max_radius: params.radius,
        };
        for attempt in 0..100000000 {
            if let Some(newe) = lj.plan_move(&mut rng, 0.03 * units::SIGMA) {
                if newe < lj.E {
                    println!("reduced energy to {:?} after {} attempts", newe, attempt);
                    lj.confirm();
                }
                if lj.E < 0.0 * units::EPSILON {
                    return lj;
                }
            }
        }
        println!("after some relaxing energy {}", lj.E);
        lj
    }
}
impl From<GrandLjParams> for Lj {
    fn from(params: GrandLjParams) -> Lj {
        Lj {
            E: 0.0 * units::EPSILON,
            error: 0.0 * units::EPSILON,
            possible_change: Change::None,
            positions: Vec::new(),
            max_radius_squared: params.radius * params.radius,
            max_radius: params.radius,
        }
    }
}

impl Lj {
    fn compute_one_atom_energy(&self, i: usize) -> Energy {
        let r2 = self.positions[i];
        let mut e: Energy = units::EPSILON * 0.0;
        for (j, &r1) in self.positions.iter().enumerate() {
            if i != j {
                e += potential((r1 - r2).norm2());
            }
        }
        e
    }
}

impl System for Lj {
    fn energy(&self) -> Energy {
        self.E
    }
    fn compute_energy(&self) -> Energy {
        let mut e: Energy = units::EPSILON * 0.0;
        for (which, &r1) in self.positions.iter().enumerate() {
            for r2 in self.positions.iter().take(which).cloned() {
                e += potential((r1 - r2).norm2());
            }
        }
        e
    }
    fn lowest_possible_energy(&self) -> Option<Energy> {
        let n = self.positions.len() as f64;
        Some(-0.5 * n * (n - 1.0) * units::EPSILON)
    }
    fn verify_energy(&self) {
        let egood = self.compute_energy();
        let expected = self.expected_accuracy(self.E);
        if (egood - self.E).abs() > expected {
            println!(
                "Error in E is {} when it should be {} < {}",
                PrettyFloat(*((egood - self.E) / units::EPSILON).value()),
                PrettyFloat(*(self.error / units::EPSILON).value()),
                PrettyFloat(*(expected / units::EPSILON).value())
            );
            assert_eq!(egood, self.E);
        }
    }
    fn randomize(&mut self, rng: &mut MyRng) -> Energy {
        for x in self.positions.iter_mut() {
            let mut r;
            loop {
                r = Vector3d::new(
                    rng.sample(Uniform::new(-1.0, 1.0)),
                    rng.sample(Uniform::new(-1.0, 1.0)),
                    rng.sample(Uniform::new(-1.0, 1.0)),
                );
                if r.norm2() < 1.0 {
                    break;
                }
            }
            *x = r*self.max_radius;
        }
        self.E = self.compute_energy();
        self.E
    }
}

impl GrandSystem for Lj {
    fn plan_add(&mut self, rng: &mut MyRng) -> Option<Energy> {
        let mut r;
        loop {
            r = Vector3d::new(
                rng.sample(Uniform::new(-1.0, 1.0)),
                rng.sample(Uniform::new(-1.0, 1.0)),
                rng.sample(Uniform::new(-1.0, 1.0)),
            );
            if r.norm2() < 1.0 {
                break;
            }
        }
        self.positions.push(r * self.max_radius);
        let e = self.E + self.compute_one_atom_energy(self.positions.len() - 1);
        let to = self.positions.pop().unwrap();
        self.possible_change = Change::Add { to, e };
        Some(e)
    }
    fn plan_remove(&mut self, rng: &mut MyRng) -> Energy {
        let which = rng.sample(Uniform::new(0, self.positions.len()));
        let e = self.E - self.compute_one_atom_energy(which);
        self.possible_change = Change::Remove { which, e };
        e
    }
    fn num_atoms(&self) -> usize {
        self.positions.len()
    }
}

impl GrandReplicaSystem for Lj {
    fn plan_swap_atom(&self, other: &Self, rng: &mut MyRng) -> Option<(usize, Energy, Energy)> {
        let which = rng.sample(Uniform::new(0, self.num_atoms()));
        let e_self = self.E - self.compute_one_atom_energy(which);
        let mut other: Lj = other.clone();
        other.positions.push(self.positions[which]);
        let e_other = other.E + other.compute_one_atom_energy(self.positions.len() - 1);
        Some((which, e_self, e_other))
    }
    fn swap_atom(&mut self, other: &mut Self, which: usize) {
        let r = self.positions.swap_remove(which);
        self.E = self.compute_energy();
        other.positions.push(r);
        other.E = other.compute_energy();
    }
}

impl ConfirmSystem for Lj {
    fn confirm(&mut self) {
        match self.possible_change {
            Change::None => (),
            Change::Move { which, to, e } => {
                self.positions[which] = to;
                self.possible_change = Change::None;
                self.set_energy(e);
            }
            Change::Add { to, e } => {
                self.positions.push(to);
                self.possible_change = Change::None;
                self.set_energy(e);
            }
            Change::Remove { which, e } => {
                self.positions.swap_remove(which);
                self.possible_change = Change::None;
                self.set_energy(e);
            }
        }
    }
    fn describe(&self) -> String {
        format!("N = {} ", self.num_atoms())
    }
}

impl MovableSystem for Lj {
    fn plan_move(&mut self, rng: &mut MyRng, mean_distance: Length) -> Option<Energy> {
        use crate::rng::vector;
        if self.positions.len() > 0 {
            let which = rng.sample(Uniform::new(0, self.positions.len()));
            let to = unsafe { *self.positions.get_unchecked(which) } + vector(rng) * mean_distance;
            self.move_atom(which, to)
        } else {
            None
        }
    }
}

#[test]
fn energy_is_right_n3() {
    energy_is_right(3);
}
#[test]
fn energy_is_right_n50() {
    energy_is_right(50);
}
#[test]
fn energy_is_right_n100() {
    energy_is_right(100);
}
#[test]
fn energy_is_right_n200() {
    energy_is_right(200);
}

#[cfg(test)]
fn energy_is_right(natoms: usize) {
    let mut lj = mk_lj(natoms);
    assert_eq!(lj.energy(), lj.compute_energy());
    let mut rng = MyRng::seed_from_u64(1);
    let mut old_energy = lj.energy();
    let maxe = (natoms as f64) * 16.0 * units::EPSILON;
    let mut i = 0.0;
    while i < 1000.0 {
        if let Some(newe) = lj.plan_move(&mut rng, Length::new(1.0)) {
            if newe < maxe || newe < old_energy {
                lj.confirm();
                println!(
                    "after move {}... {} vs {}",
                    i,
                    lj.energy(),
                    lj.compute_energy()
                );
                lj.verify_energy();
                old_energy = newe;
                i += 1.0;
            } else {
                println!(
                    "rejected move giving {} (vs old_energy {} and maxe {})",
                    newe, old_energy, maxe
                );
                i += 1e-6;
            }
        }
    }
}

#[cfg(test)]
fn mk_lj(natoms: usize) -> Lj {
    let radius = 2.0 * (natoms as f64).powf(1.0 / 3.0) * units::SIGMA;
    let radius = 5.0 * radius;
    Lj::from(LjParams {
        N: natoms,
        radius,
    })
}

#[test]
fn init_lj() {
    mk_lj(50);
}

#[test]
fn test_swap() {
    let mut system_one = mk_lj(5);
    let mut system_two = mk_lj(6);
    let mut ran = MyRng::seed_from_u64(1);
    let (which, energy_two, energy_one) = system_two.plan_swap_atom(&system_one, &mut ran).unwrap();
    system_two.swap_atom(&mut system_one, which);
    assert!((system_one.energy() - energy_one).abs() < Energy::new(0.00001));
    assert!((system_two.energy() - energy_two).abs() < Energy::new(0.00001));
    assert_eq!(system_one.num_atoms(), 6);
    assert_eq!(system_two.num_atoms(), 5);
}
