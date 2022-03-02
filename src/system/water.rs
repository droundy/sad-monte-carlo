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

// alpha = 5.213 kJ/mol => polarizability const

/// Angle of water molecule
const ANGLE_HOH: f64 = 1.911; // units of radians
/// Charge of oxygen
const CHARGE_O: f64 = -0.8476; // units of elementary charges
/// Charge of hydrogen
const CHARGE_H: f64 = -CHARGE_O / 2.0; // units of elementary charges

/// The parameters needed to configure a water system.
#[allow(non_snake_case)]
#[derive(Serialize, Deserialize, Debug, AutoArgs)] // some serde -> serialize/deserialize 
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


//For potential below, induction energy or polarization energy needs to be added.
// dipole moment is charge times distance of seperation. 1 D = 3.33564e-30 (Cm) (dim=>QL)

//Ewald?

    //adding induction energy
    //let mut polarization: Energy = 0.0 * units::EPSILON;
    // (1/2a)* sum( for all j) [(mu(j)-mu0)**2]
    // a = alpha = polarizability of a water molecule
    // mu0 = dipole moment of isolated molecule = 1.85 D
    // Does this need to go outside the function with the caller?
    // sum over what?

    //let mut i = 0;
   // sum thing


   // FIX BAD units 
//const ALPHA: f64 = 5.213; // kJ/mol ..alpha => polirizability 
//const MU: f64 = 0.415; //2.3-1.885 = 1.115. Units are D, 1D = 3.335 xx 10^-30 Coulumbs * meters
                        // 2.3 D - 1.885 D= 0.415 = mu-liquid - mu-gas => liquid and gas phase dipole moment
                        // 2.43 D to 3.09 D have been used supposedly.


/// Compute the energy between two molecules. SPC/E
fn potential(a: WaterMolecule, b: WaterMolecule) -> Energy { // Returns an energy *
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

    let sigma = 3.1656; // Angstroms => L
    //let epsilon_LJ = 78.24; // * KBoltzman  
    //let k_boltzman = 1.3806e-23; // J/K => Q/T
    let epsilon_lj_k_boltzman = 1.0802e-21; // J/K => Energy/Temp

    // FIX ME bad units on epsilon_lj_k_boltzman ?
    let sigma_over_r_2 = units::SIGMA * units::SIGMA *sigma*sigma / (a.position - b.position).norm2(); // sigma => length
    let sigma_over_r_6 = sigma_over_r_2 * sigma_over_r_2 * sigma_over_r_2;
    let sigma_over_r_12 = sigma_over_r_6 * sigma_over_r_6;
    potential += 4.0 * units::EPSILON * epsilon_lj_k_boltzman * (sigma_over_r_12 - sigma_over_r_6); // 4 * epsilon LJ = 4*78.24*(1.38e-23 J/molK) =>4*78.24*Kb *(....12-6)?

    // Does this make sense here? NO Outside with energy callers.. 
    //move atom shouldn't call this.
    //potential += units::EPSILON*(MU*MU) / ( 2.0 * ALPHA);
    potential
}







/// Compute electric potential energy between two point charges.
fn electric_potential(a: Vector3d<Length>, b: Vector3d<Length>, qq: f64) -> Energy {
    let coulomb_constant = 675.109691 * units::EPSILON * units::SIGMA; // x = ky => x/k = y => y = 675.109691(Newton??? * Angstrom^2/elementary charge^2)/(9e9(N*m*m/C*C)) = 7.50121879e-8()
    return coulomb_constant * qq / (a - b).norm2().sqrt(); // What are the units here? meV?
}                   // qq has no units built in, but is elementary charge^2          

// E => Nm =J... Kgm^2/s^2
// is the mass of a proton the unit of mass
// mp = 1.67e-27 kg
// E => mp*A^2/s^2 ? 
// kqq/r => Nmm/CC
// Kgm^3/(Cs)2 *(AMU/Kg) *(A/m)^3 *(C/e)^2
// e=> elementary charge, A => angstrom


impl Water {
    /// Change the position of a specified particle.  Returns the change in energy, or
    /// `None` if the particle could not be placed there.
    pub fn move_atom(&mut self, which: usize, r: Vector3d<Length>) -> Option<Energy> { // Why is there an Option? It could either return an energy, or none.
        let old = self.molecules[which]; //                                             which is just a number, type usize "The size of this primitive is how many bytes it takes to reference any location in memory."
        let previous_rsqr = old.position.norm2();
        if r.norm2() > self.max_radius_squared && r.norm2() > previous_rsqr { // if r^2 is outside The square of the maximum radius permitted
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
        let zero_molecule = WaterMolecule {
            position: zero_vector,
            h1: zero_vector,
            h2: zero_vector,
        };
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

#[allow(non_snake_case)]
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

        //let mut mu: f64 = 0.415; // D, 1D = 3.335e-30 Couloumbs * meters... 2.3 D - 1.885 D= 0.415
        // Debyes to Couloumb meters to elementary charge * angstroms
        //
        //let couloumbMeter_debye = 3.335e-30; 
        //let elementaryCharge_Couloumb = 1.6022e-19;
        //let angstrom_meter = 1.0e10;
        //mu = mu * couloumbMeter_debye *elementaryCharge_Couloumb *angstrom_meter;
        //Computed the above for mu in python:
        //FIX ME bad units on alpha... mu (fixed)
        let alpha: f64 = 5213.0; // polirizability => J/mol 
        let mu = 2.2175e-39; //elementary Charge * angstroms
        let mMolar_H20: f64 = 18.015; // moles (15.999 moles O + 2*1.008 moles H)
        e += units::EPSILON*(mu*mu) / ( 2.0 * alpha * mMolar_H20 *self.molecules.len() as f64 );// FIX ME bad units still
        e

        //alpha units => J/mol = Nm/mol what is the energy unit we get from potential
    }

    //Can't add... undeclared func for system in mod.rs
    //don't even need to call elsewhere as this is specific to this energy
    // fn polarizationEnergy(&self) -> Energy {
    //     let mut ePol: Energy = units::EPSILON * 0.0;
    //     ePol
    // }

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
    fn dimensionality(&self) -> u64 {
        self.min_moves_to_randomize() * 3
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
            let to = unsafe { *self.molecules.get_unchecked(which) }.position
                + vector(rng) * mean_distance;
            self.move_atom(which, to)
        } else {
            None
        }
    }
    fn max_size(&self) -> Length {
        self.max_radius
    }
}

/// Generate a random number uniformly in `[a, b)`. FIXME eliminate
fn rand_uniform(rng: &mut MyRng, a: f64, b: f64) -> f64 {
    rng.sample(Uniform::new(a, b))
}

/// Generate a random vector uniformly in the unit ball.
pub fn rand_unit_ball(rng: &mut MyRng) -> Vector3d<f64> {
    loop {
        let r = Vector3d::new(
            //rng.gen_range(-1.0, 1.0),
            //rng.gen_range(-1.0, 1.0),
            rand_uniform(rng, -1.0, 1.0),
            rand_uniform(rng, -1.0, 1.0),
            rand_uniform(rng, -1.0, 1.0),
        );
        if r.norm2() < 1.0 {
            return r;
        }
    }
}

/// Generate a random vector uniformly in the unit ball.
pub fn rand_unit_ball_2(rng: &mut MyRng) -> Vector3d<f64> {
    loop {
        let r = Vector3d::new(
            rng.gen_range(-1.0, 1.0),
            rng.gen_range(-1.0, 1.0),
            rng.gen_range(-1.0, 1.0),
        );
        if r.norm2() < 1.0 {
            return r;
        }
    }
}

/// Generate a random 3D rotation.  FIXME move to the rotation module as a
/// constructor method
fn rand_rotation(rng: &mut MyRng) -> Rotation {
    let (x1, y1) = rand_uniform(rng, 0.0, 2.0 * std::f64::consts::PI).sin_cos();
    let (x2, y2) = rand_uniform(rng, 0.0, 2.0 * std::f64::consts::PI).sin_cos();
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
    let maxe = (natoms as f64) * 16.0 * units::EPSILON; // Why 16 * natoms?
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
    Water::from(WaterParams { N: natoms, radius }) // Sets params
}

#[test]
fn init_water() {
    mk_water(50);
}

