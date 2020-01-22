//! A square well fluid.

use super::*;

use dimensioned::{Dimensionless, Abs, Sqrt};
use vector3d::Vector3d;
use rand::prelude::*;
use rand::distributions::Uniform;

/// The parameters needed to configure a Weeks-Chandler-Anderson (WCA) system.
#[allow(non_snake_case)]
#[derive(Serialize, Deserialize, Debug, ClapMe)]
pub struct LjParams {
    /// The number of atoms
    N: usize,
    /// The radius of the spherical box
    radius: Length,
    /// The number of bins for the radial distribution
    n_radial: Option<usize>,
}

/// The parameters needed to configure a Weeks-Chandler-Anderson (WCA) system for a grand computation.
#[allow(non_snake_case)]
#[derive(Serialize, Deserialize, Debug, ClapMe)]
pub struct GrandLjParams {
    /// The radius of the spherical box
    radius: Length,
    /// The number of bins for the radial distribution
    n_radial: Option<usize>,
}

#[allow(non_snake_case)]
/// A WCA fluid.
#[derive(Serialize, Deserialize, Debug)]
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
    /// Number of points in radial distribution
    n_radial: Option<usize>,
}

/// This defines the energy/radial bins.
#[derive(Serialize, Deserialize, Debug, Default)]
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
    Move { which: usize, to: Vector3d<Length>, e: Energy },
    Add { to: Vector3d<Length>, e: Energy },
    Remove { which: usize, e: Energy },
    /// Make no changes to the system
    None,
}

/// Define the WCA interaction potential and criteria
fn potential(r_squared: Area) -> Energy {
    let sig_sqr = units::SIGMA*units::SIGMA;
    4.0*units::EPSILON*((sig_sqr/r_squared).powi(6) - (sig_sqr/r_squared).powi(3))
}

impl Lj {
    fn cm(&self) -> Vector3d<Length> {
        let mut cm = self.positions[0];
        for x in self.positions.iter().skip(1) {
            cm = cm + *x;
        }
        cm/self.positions.len() as f64
    }
    /// Move a specified atom.  Returns the change in energy, or
    /// `None` if the atom could not be placed there.
    pub fn move_atom(&mut self, which: usize, r: Vector3d<Length>) -> Option<Energy> {
        let previous_rsqr = self.positions[which].norm2();
        if r.norm2() > self.max_radius_squared && r.norm2() > previous_rsqr {
            return None;
        }
        let mut e = self.E;
        let from = self.positions[which];
        for r1 in self.positions.iter().cloned().enumerate().filter(|&(i,_)| i != which).map(|(_,x)| x) {
            e += potential((r1-r).norm2()) - potential((r1-from).norm2());
        }
        self.possible_change = Change::Move{ which, to: r, e };
        Some(e)
    }
    fn expected_accuracy(&self, newe: Energy) -> Energy {
        newe.abs()*1e-14*(self.positions.len() as f64)*(self.positions.len() as f64)
    }

    fn set_energy(&mut self, new_e: Energy) {
        let new_error = if new_e.abs() > self.E.abs() {
            new_e.abs()*1e-15*(self.positions.len() as f64)
        } else {
            self.E.abs()*1e-15*(self.positions.len() as f64)
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
        let mut rng = crate::rng::MyRng::from_u64(0);
        let mut best_energy = 1e80*units::EPSILON;
        let mut best_positions = Vec::new();
        println!("I am creating a LJ system with {} atoms!", params.N);
        for attempt in 0..10000000 {
            let mut positions = Vec::new();
            for _ in 0..params.N {
                let mut r;
                loop {
                    r = Vector3d::new(rng.sample(Uniform::new(-1.0, 1.0)),
                                      rng.sample(Uniform::new(-1.0, 1.0)),
                                      rng.sample(Uniform::new(-1.0, 1.0)),);
                    if r.norm2() < 1.0 {
                        break;
                    }
                }
                positions.push(r*params.radius);
            }
            let mut cm = positions[0];
            for x in positions.iter().skip(1) {
                cm = cm + *x;
            }
            cm = cm/params.N as f64;
            for x in positions.iter_mut() {
                *x = *x - cm;
            }
            if positions.iter().any(|x| x.norm2() > params.radius*params.radius) {
                continue;
            }
            let mut lj = Lj {
                E: 0.0*units::EPSILON,
                error: 0.0*units::EPSILON,
                possible_change: Change::None,
                positions,
                max_radius_squared: params.radius*params.radius,
                max_radius: params.radius,
                n_radial: params.n_radial,
            };
            lj.E = lj.compute_energy();
            if lj.E < best_energy {
                best_energy = lj.E;
                best_positions = lj.positions.clone();
                println!("found a new best energy{:?} after {} attempts", lj.E, attempt);
            }
            if lj.E < 0.0*units::EPSILON {
                return lj;
            }
        }
        let mut lj = Lj {
            E: best_energy,
            error: 0.0*units::EPSILON,
            possible_change: Change::None,
            positions: best_positions,
            max_radius_squared: params.radius*params.radius,
            max_radius: params.radius,
            n_radial: params.n_radial,
        };
        for attempt in 0..100000000 {
            if let Some(newe) = lj.plan_move(&mut rng, 0.03*units::SIGMA) {
                if newe < lj.E {
                    println!("reduced energy to {:?} after {} attempts", newe, attempt);
                    lj.confirm();
                }
                if lj.E < 0.0*units::EPSILON {
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
            E: 0.0*units::EPSILON,
            error: 0.0*units::EPSILON,
            possible_change: Change::None,
            positions: Vec::new(),
            max_radius_squared: params.radius*params.radius,
            max_radius: params.radius,
            n_radial: params.n_radial,
        }
    }
}

impl Lj {
    fn compute_one_atom_energy(&self, i: usize) -> Energy {
        let r2 = self.positions[i];
        let mut e: Energy = units::EPSILON*0.0;
        for (j, &r1) in self.positions.iter().enumerate() {
            if i != j {
                e += potential((r1-r2).norm2());
            }
        }
        e
    }
}

impl System for Lj {
    type CollectedData = Collected;
    fn collect_data(&self, data: &mut Collected, iter: u64) {
        if iter % ((self.positions.len()*self.positions.len()) as u64) == 0 {
            if let Some(n) = self.n_radial {
                if iter as usize  % self.positions.len() == 0 {
                    if data.from_center.len() != n {
                        data.from_center = vec![0; n];
                    }
                    if data.from_cm.len() != n {
                        data.from_cm = vec![0; n];
                    }
                    let cm = self.cm();
                    let n = n as f64;
                    for x in self.positions.iter().cloned() {
                        let r = x.norm2().sqrt();
                        let i_r = (*(r/self.max_radius).value() * n).floor() as usize;
                        data.from_center[i_r] += 1;

                        let r = (x-cm).norm2().sqrt();
                        if r < self.max_radius {
                            let i_r = (*(r/self.max_radius).value() * n).floor() as usize;
                            data.from_cm[i_r] += 1;
                        }
                    }
                }
            }
        }
    }
    fn energy(&self) -> Energy {
        self.E
    }
    fn compute_energy(&self) -> Energy {
        let mut e: Energy = units::EPSILON*0.0;
        for (which, &r1) in self.positions.iter().enumerate() {
            for r2 in self.positions.iter().take(which).cloned() {
                e += potential((r1-r2).norm2());
            }
        }
        e
    }
    fn update_caches(&mut self) {
    }
    fn lowest_possible_energy(&self) -> Option<Energy> {
        let n = self.positions.len() as f64;
        Some(-0.5*n*(n-1.0)*units::EPSILON)
    }
    fn verify_energy(&self) {
        let egood = self.compute_energy();
        let expected = self.expected_accuracy(self.E);
        if (egood - self.E).abs() > expected {
            println!("Error in E is {} when it should be {} < {}", egood - self.E, self.error, expected);
            assert_eq!(egood, self.E);
        }
    }
}

impl GrandSystem for Lj {
    fn plan_add(&mut self, rng: &mut MyRng) -> Option<Energy> {
        let mut r;
        loop {
            r = Vector3d::new(rng.sample(Uniform::new(-1.0, 1.0)),
                              rng.sample(Uniform::new(-1.0, 1.0)),
                              rng.sample(Uniform::new(-1.0, 1.0)),);
            if r.norm2() < 1.0 {
                break;
            }
        }
        self.positions.push(r*self.max_radius);
        let e = self.E + self.compute_one_atom_energy(self.positions.len()-1);
        let to = self.positions.pop().unwrap();
        self.possible_change = Change::Add { to, e, };
        Some(e)
    }
    fn plan_remove(&mut self, rng: &mut MyRng) -> Energy {
        let which = rng.sample(Uniform::new(0, self.positions.len()));
        let e = self.E - self.compute_one_atom_energy(which);
        self.possible_change = Change::Remove { which, e, };
        e
    }
    fn num_atoms(&self) -> usize {
        self.positions.len()
    }
}

impl ConfirmSystem for Lj {
    fn confirm(&mut self) {
        match self.possible_change {
            Change::None => (),
            Change::Move{which, to, e} => {
                self.positions[which] = to;
                self.possible_change = Change::None;
                self.set_energy(e);
            },
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
}

impl MovableSystem for Lj {
    fn plan_move(&mut self, rng: &mut MyRng, mean_distance: Length) -> Option<Energy> {
        if self.positions.len() > 0 {
            let which = rng.sample(Uniform::new(0, self.positions.len()));
            let to = unsafe { *self.positions.get_unchecked(which) } + rng.vector()*mean_distance;
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
    let mut rng = MyRng::from_u64(1);
    let mut old_energy = lj.energy();
    let maxe = (natoms as f64)*16.0*units::EPSILON;
    let mut i = 0.0;
    while i < 1000.0 {
        if let Some(newe) = lj.plan_move(&mut rng, Length::new(1.0)) {
            if newe < maxe || newe < old_energy {
                lj.confirm();
                println!("after move {}... {} vs {}", i, lj.energy(), lj.compute_energy());
                lj.verify_energy();
                old_energy = newe;
                i += 1.0;
            } else {
                println!("rejected move giving {} (vs old_energy {} and maxe {})", newe, old_energy, maxe);
                i += 1e-6;
            }
        }
    }
}

#[cfg(test)]
fn mk_lj(natoms: usize) -> Lj {
    let radius = 2.0*(natoms as f64).powf(1.0/3.0)*units::SIGMA;
    let radius = 5.0*radius;
    Lj::from(LjParams { N: natoms, radius, n_radial: None })
}

#[test]
fn init_lj() {
    mk_lj(50);
}
