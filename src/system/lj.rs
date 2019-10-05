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
    /// Energy bin width for radial info
    width: Energy,
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
    /// Histogram of radial distribution
    bins: Bins,
}

/// This defines the energy/radial bins.
#[derive(Serialize, Deserialize, Debug)]
pub struct Bins {
    /// The lowest allowed energy in any bin.
    pub min: Energy,
    /// The energy bin size.
    pub width: Energy,
    /// The ln weight for each energy bin.
    pub radial: Vec<Vec<u64>>,
    /// Number of points in radial distribution
    n_radial: usize,
}

impl Bins {
    /// Find the index corresponding to a given energy.  This should
    /// panic if the energy is less than `min`.
    pub fn energy_to_index(&self, s: Energy) -> usize {
        *((s - self.min)/self.width).value() as usize
    }
    /// Find the energy corresponding to a given index.
    pub fn index_to_energy(&self, i: usize) -> Energy {
        self.min + (i as f64 + 0.5)*self.width
    }
    /// Make room in our arrays for a new energy value
    pub fn prepare_for_energy(&mut self, e: Energy) {
        assert!(self.width > Energy::new(0.0));
        while e < self.min {
            // this is a little wasteful, but seems the easiest way to
            // ensure we end up with enough room.
            self.radial.insert(0, vec![0u64; self.n_radial]);
        }
        while e >= self.min + self.width*(self.radial.len() as f64) {
            self.radial.push(vec![0u64; self.n_radial]);
        }
    }
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
/// Define the types of changes that can be made to the system
enum Change {
    /// Move an atom already in the system
    Move { which: usize, to: Vector3d<Length>, e: Energy },
    /// Make no changes to the system
    None,
}

/// Define the WCA interaction potential and criteria
fn potential(r_squared: Area) -> Energy {
    let sig_sqr = units::SIGMA*units::SIGMA;
    4.0*units::EPSILON*((sig_sqr/r_squared).powi(6) - (sig_sqr/r_squared).powi(3))
}

impl Lj {
    fn increment_radial(&mut self, r: Length) {
        let i_r = (*(r/self.max_radius).value() * self.bins.n_radial as f64).floor() as usize;
        let i_e = self.bins.energy_to_index(self.E);
        self.bins.radial[i_e][i_r] += 1;
    }
    /// Move a specified atom.  Returns the change in energy, or
    /// `None` if the atom could not be placed there.
    pub fn move_atom(&mut self, which: usize, r: Vector3d<Length>) -> Option<Energy> {
        let dr = r - self.positions[which];
        let r_from_cm = r - dr/(self.positions.len() as f64);
        let previous_rsqr = self.positions[which].norm2();
        self.increment_radial(previous_rsqr.sqrt());
        if r_from_cm.norm2() > self.max_radius_squared && r_from_cm.norm2() > previous_rsqr {
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
        let mut positions = Vec::new();

        let mut rng = ::rng::MyRng::from_u64(0);
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
        let mut lj = Lj {
            E: 0.0*units::EPSILON,
            error: 0.0*units::EPSILON,
            possible_change: Change::None,
            positions,
            max_radius_squared: params.radius*params.radius,
            max_radius: params.radius,
            bins: Bins {
                n_radial: params.n_radial.unwrap_or(100),
                radial: Vec::new(),
                min: 0.0*units::EPSILON,
                width: params.width,
            },
        };
        lj.E = lj.compute_energy();
        lj.bins.min = ((lj.E/params.width).value().round() - 0.5)*params.width;
        lj
    }
}

impl System for Lj {
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

impl ConfirmSystem for Lj {
    fn confirm(&mut self) {
        match self.possible_change {
            Change::None => (),
            Change::Move{which, to, e} => {
                self.bins.prepare_for_energy(e);
                self.positions[which] = to;
                let mut cm = self.positions[0];
                for x in self.positions.iter().skip(1) {
                    cm = cm + *x;
                }
                cm = cm/self.positions.len() as f64;
                for x in self.positions.iter_mut() {
                    *x = *x - cm;
                }
                self.possible_change = Change::None;
                self.set_energy(e);
            },
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
    Lj::from(LjParams { N: natoms, radius })
}

#[test]
fn init_lj() {
    mk_lj(50);
}
