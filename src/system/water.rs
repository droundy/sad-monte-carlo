//! A bunch of water.

use super::*;

use crate::prettyfloat::PrettyFloat;
use crate::rotation::Rotation;
use crate::unit_quaternion::UnitQuaternion;
use dimensioned::{Abs, Dimensionless, Sqrt};
use rand::distributions::Uniform;
use rand::{Rng, SeedableRng};
use vector3d::Vector3d;

// sigma = 3.166 Å
// epsilon = 6.737 meV

/// Angle of water molecule
const ANGLE_HOH: f64 = 1.911; // units of radians
/// Charge of oxygen
const CHARGE_O: f64 = -0.8476; // units of elementary charges
/// Charge of hydrogen
const CHARGE_H: f64 = -CHARGE_O / 2.0; // units of elementary charges

/// The parameters needed to configure a water system.
#[allow(non_snake_case)]
#[derive(Serialize, Deserialize, Debug, AutoArgs)]
pub struct WaterParams {
    /// The number of molecules
    N: usize,
    /// The radius of the spherical box
    radius: Length,
}

/// Some water.
#[allow(non_snake_case)]
#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct Water {
    /// The energy of the system
    E: Energy,
    /// The estimated accumulated error so far in E
    error: Energy,
    /// The last change we made (and might want to undo).
    possible_change: Change,
    /// The molecules
    molecules: Vec<WaterMolecule>,
    /// The square of the maximum radius permitted
    max_radius_squared: Area,
    /// The maximum radius permitted
    max_radius: Length,
}

/// A single water.
#[derive(Copy, Clone, Serialize, Deserialize, Debug)]
struct WaterMolecule {
    // Maybe Vec<Vector3d<Length>>?  :(
    /// The position of the oxygen.
    position: Vector3d<Length>,
    /// The vector from the oxygen to the first hydrogen.
    h1: Vector3d<Length>,
    /// The vector from the oxygen to the second hydrogen.
    h2: Vector3d<Length>,
}

/// This defines the energy/radial bins.
#[derive(Serialize, Deserialize, Debug, Default, Clone)]
pub struct Collected {
    /// The number of atoms for each radial bin from center of sphere.
    pub from_center: Vec<u64>,
    /// The number of atoms for each radial bin from center of mass.
    pub from_cm: Vec<u64>,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
/// Define the types of changes that can be made to the system
enum Change {
    /// Move and rotate a particle already in the system
    Move {
        which: usize,
        new: WaterMolecule,
        e: Energy,
    },
    /// Make no changes to the system
    None,
}

/// Compute the energy between two molecules. SPC/E
fn potential(a: WaterMolecule, b: WaterMolecule) -> Energy {
    let mut potential: Energy = 0.0 * units::EPSILON;
    potential += electric_potential(a.position, b.position, CHARGE_O * CHARGE_O);
    potential += electric_potential(a.position, b.h1, CHARGE_O * CHARGE_H);
    potential += electric_potential(a.position, b.h2, CHARGE_O * CHARGE_H);
    potential += electric_potential(a.h1, b.position, CHARGE_H * CHARGE_O);
    potential += electric_potential(a.h1, b.h1, CHARGE_H * CHARGE_H);
    potential += electric_potential(a.h1, b.h2, CHARGE_H * CHARGE_H);
    potential += electric_potential(a.h2, b.position, CHARGE_H * CHARGE_O);
    potential += electric_potential(a.h2, b.h1, CHARGE_H * CHARGE_H);
    potential += electric_potential(a.h2, b.h2, CHARGE_H * CHARGE_H);
    let sigma_over_r_2 = units::SIGMA * units::SIGMA / (a.position - b.position).norm2();
    let sigma_over_r_6 = sigma_over_r_2 * sigma_over_r_2 * sigma_over_r_2;
    let sigma_over_r_12 = sigma_over_r_6 * sigma_over_r_6;
    potential += 4.0 * units::EPSILON * (sigma_over_r_12 - sigma_over_r_6);
    potential
}

/// Compute electric potential energy between two point charges.
fn electric_potential(a: Vector3d<Length>, b: Vector3d<Length>, qq: f64) -> Energy {
    let coulomb_constant = 675.109691 * units::EPSILON * units::SIGMA;
    return coulomb_constant * qq / (a - b).norm2().sqrt();
}

impl Water {
    /// Change the position of a specified particle.  Returns the change in energy, or
    /// `None` if the particle could not be placed there.
    pub fn move_atom(&mut self, which: usize, r: Vector3d<Length>) -> Option<Energy> {
        let old = self.molecules[which];
        let previous_rsqr = old.position.norm2();
        if r.norm2() > self.max_radius_squared && r.norm2() > previous_rsqr {
            return None;
        }
        let mut e = self.E;
        let new = WaterMolecule { position: r, ..old };
        for molecule in self
            .molecules
            .iter()
            .cloned()
            .enumerate()
            .filter(|&(i, _)| i != which)
            .map(|(_, x)| x)
        {
            e += potential(molecule, new) - potential(molecule, old);
        }
        self.possible_change = Change::Move { which, new, e };
        Some(e)
    }
    fn expected_accuracy(&self, newe: Energy) -> Energy {
        // TODO: where does this formula come from
        newe.abs() * 1e-14 * (self.molecules.len() as f64) * (self.molecules.len() as f64)
    }
    fn set_energy(&mut self, new_e: Energy) {
        // TODO: what is this
        let new_error = if new_e.abs() > self.E.abs() {
            new_e.abs() * 1e-15 * (self.molecules.len() as f64)
        } else {
            self.E.abs() * 1e-15 * (self.molecules.len() as f64)
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

impl From<WaterParams> for Water {
    fn from(params: WaterParams) -> Water {
        let mut rng = crate::rng::MyRng::seed_from_u64(0);
        let zero_vector = Vector3d::new(0.0, 0.0, 0.0) * units::SIGMA;
        let zero_molecule = WaterMolecule { position: zero_vector, h1: zero_vector, h2: zero_vector };
        let mut water = Water {
            E: 0.0 * units::EPSILON,
            error: 0.0 * units::EPSILON,
            possible_change: Change::None,
            molecules: vec![zero_molecule; params.N],
            max_radius_squared: params.radius * params.radius,
            max_radius: params.radius,
        };
        water.randomize(&mut rng);
        water
    }
}

impl System for Water {
    fn energy(&self) -> Energy {
        self.E
    }
    fn compute_energy(&self) -> Energy {
        let mut e: Energy = units::EPSILON * 0.0;
        for (which, &m1) in self.molecules.iter().enumerate() {
            for m2 in self.molecules.iter().take(which).cloned() {
                e += potential(m1, m2);
            }
        }
        e
    }
    // fn lowest_possible_energy(&self) -> Option<Energy> {
    //     // TODO: where does this formula come from
    //     let n = self.molecules.len() as f64;
    //     Some(-0.5 * n * (n - 1.0) * units::EPSILON)
    // }
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
        for x in self.molecules.iter_mut() {
            *x = rand_molecule(rng, self.max_radius);
        }
        self.E = self.compute_energy();
        self.E
    }
    fn min_moves_to_randomize(&self) -> u64 {
        self.molecules.len() as u64
    }
}

impl ConfirmSystem for Water {
    fn confirm(&mut self) {
        match self.possible_change {
            Change::None => (),
            Change::Move { which, new, e } => {
                self.molecules[which] = new;
                self.possible_change = Change::None;
                self.set_energy(e);
            }
        }
    }
    fn describe(&self) -> String {
        format!("N = {} ", self.molecules.len())
    }
}

impl MovableSystem for Water {
    fn plan_move(&mut self, rng: &mut MyRng, mean_distance: Length) -> Option<Energy> {
        use crate::rng::vector;
        if self.molecules.len() > 0 {
            let which = rng.sample(Uniform::new(0, self.molecules.len()));
            let to = unsafe { *self.molecules.get_unchecked(which) }.position + vector(rng) * mean_distance;
            self.move_atom(which, to)
        } else {
            None
        }
    }
}

/// Generate a random number uniformly in `[a, b)`. FIXME eliminate
fn rand_uniform(rng: &mut MyRng, a: f64, b: f64) -> f64 {
    rng.sample(Uniform::new(a, b))
}

/// Generate a random vector uniformly in the unit ball.
fn rand_unit_ball(rng: &mut MyRng) -> Vector3d<f64> {
    loop {
        let r = Vector3d::new(
            rng.gen_range(-1.0, 1.0),
            rand_uniform(rng, -1.0, 1.0),
            rand_uniform(rng, -1.0, 1.0),
        );
        if r.norm2() < 1.0 {
            return r;
        }
    }
}

/// Generate a random 3D rotation.  FIXME move to the rotation module as a
/// constructor method
fn rand_rotation(rng: &mut MyRng) -> Rotation {
    let (x1, y1) = rand_uniform(rng, 0.0, 2.0*std::f64::consts::PI).sin_cos();
    let (x2, y2) = rand_uniform(rng, 0.0, std::f64::consts::TAU).sin_cos();
    let u = rand_uniform(rng, 0.0, 1.0);
    let u1 = u.sqrt();
    let u2 = (1.0 - u).sqrt();
    let q = UnitQuaternion {
        x: u1 * x1,
        y: u1 * y1,
        z: u2 * x2,
        w: u2 * y2,
    };
    Rotation { q }
}

/// Generate a random molecule position and orientation. FIXME method WaterMolecule::random(&mut rng)
fn rand_molecule(rng: &mut MyRng, radius: Length) -> WaterMolecule {
    // Distance between oxygen and hydrogen
    let r_oh: Length = 0.315856 * units::SIGMA;
    let position = rand_unit_ball(rng) * radius;
    let rotation = rand_rotation(rng);
    let h1 = rotation * Vector3d::new(1.0, 0.0, 0.0) * r_oh;
    let h2 = rotation * Vector3d::new(ANGLE_HOH.cos(), ANGLE_HOH.sin(), 0.0) * r_oh;
    WaterMolecule { position, h1, h2 }
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
    let mut water = mk_water(natoms);
    assert_eq!(water.energy(), water.compute_energy());
    let mut rng = MyRng::seed_from_u64(1);
    let mut old_energy = water.energy();
    let maxe = (natoms as f64) * 16.0 * units::EPSILON;
    let mut i = 0.0;
    while i < 1000.0 {
        if let Some(newe) = water.plan_move(&mut rng, Length::new(1.0)) {
            if newe < maxe || newe < old_energy {
                water.confirm();
                println!(
                    "after move {}... {} vs {}",
                    i,
                    water.energy(),
                    water.compute_energy()
                );
                water.verify_energy();
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
fn mk_water(natoms: usize) -> Water {
    let radius = 2.0 * (natoms as f64).powf(1.0 / 3.0) * units::SIGMA;
    let radius = 5.0 * radius;
    Water::from(WaterParams {
        N: natoms,
        radius,
    })
}

#[test]
fn init_water() {
    mk_water(50);
}

