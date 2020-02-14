//! A Weeks-Chandler-Anderson fluid.

use super::*;

use dimensioned::{Dimensionless, Sqrt, Abs};
use vector3d::Vector3d;
use rand::prelude::*;
use rand::distributions::Uniform;
use std::default::Default;

use super::optcell::{Cell, CellDimensions};


/// The parameters needed to configure a Weeks-Chandler-Anderson (WCA) system.
#[derive(Serialize, Deserialize, Debug, ClapMe)]
pub struct WcaParams {
   _dim: CellDimensions,
}

#[allow(non_snake_case)]
/// A WCA fluid.
#[derive(Serialize, Deserialize, Debug)]
pub struct Wca {
    /// The energy of the system
    E: Energy,
    /// The estimated accumulated error so far in E
    error: Energy,
    /// The dimensions of the box.
    pub cell: Cell,
    /// The last change we made (and might want to undo).
    possible_change: Change,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
/// Define the types of changes that can be made to the system
enum Change {
    /// Move an atom already in the system
    Move { which: usize, to: Vector3d<Length>, e: Energy },
    /// Add an atom to the system
    /// The added atom is always pushed to the end of the vector!
    Add { to: Vector3d<Length>, e: Energy },
    /// Remove an atom from the system
    Remove { which: usize, e: Energy },
    /// Make no changes to the system
    None,
}

fn r_cutoff() -> Length {
    2.0_f64.powf(1.0/6.0)*units::SIGMA
}

/// Define the WCA interaction potential and criteria
fn potential(r_squared: Area) -> Energy {
    let r_cutoff: Length = r_cutoff();
    let r_cutoff_squared: Area = r_cutoff*r_cutoff;
    let sig_sqr = units::SIGMA*units::SIGMA;
    if r_squared < r_cutoff_squared {
        4.0*units::EPSILON*((sig_sqr/r_squared).powi(6) - (sig_sqr/r_squared).powi(3)) + units::EPSILON
    } else {
        0.0*units::EPSILON
    }
}

/// Define the WCA interaction pressure
fn potential_pressure(r_squared: Area) -> Energy {
    let r_cutoff: Length = r_cutoff();
    let r_cutoff_squared: Area = r_cutoff*r_cutoff;
    let sig_sqr = units::SIGMA*units::SIGMA;
    if r_squared < r_cutoff_squared {
        // This is dE/dr*r/2. the extra factor of two comes from
        // avoiding double counting.
        4.0*3.0*units::EPSILON*(-2.0*(sig_sqr/r_squared).powi(6) + (sig_sqr/r_squared).powi(3))
    } else {
        0.0*units::EPSILON
    }
}

impl Wca {
    /// Add an atom at a given location.  Returns the change in
    /// energy, or `None` if the atom could not be placed there.
    pub fn add_atom_at(&mut self, r: Vector3d<Length>) -> Option<Energy> {
        let mut e = self.E;
        for r1 in self.cell.maybe_interacting_atoms(r) {
            let dist2 = (r1-r).norm2();
            if dist2 < self.cell.r_cutoff*self.cell.r_cutoff {
                e += potential(dist2);
            }
        }
        self.possible_change = Change::Add{ to: r, e };
        Some(e)
    }
    /// Move a specified atom.  Returns the change in energy, or
    /// `None` if the atom could not be placed there.
    pub fn move_atom(&mut self, which: usize, r: Vector3d<Length>) -> Option<Energy> {
        let mut e = self.E;
        let from = self.cell.positions[which];
        for r1 in self.cell.maybe_interacting_atoms_excluding(r, which) {
            e += potential((r1-r).norm2());
        }
        for r1 in self.cell.maybe_interacting_atoms_excluding(from, which) {
            e -= potential((r1-from).norm2());
        }
        self.possible_change = Change::Move{ which, to: r, e };
        Some(e)
    }
    /// Plan to remove the specified atom.  Returns the change in energy.
    pub fn remove_atom_number(&mut self, which: usize) -> Energy {
        let r = self.cell.positions[which];
        let mut e = self.E;
        for r1 in self.cell.maybe_interacting_atoms_excluding(r, which) {
            e -= potential((r1-r).norm2());
        }
        self.possible_change = Change::Remove{ which, e };
        e
    }
    fn set_energy(&mut self, new_e: Energy) {
        let new_error = if new_e.abs() > self.E.abs() {
            new_e.abs()*1e-14*(self.num_atoms() as f64)
        } else {
            self.E.abs()*1e-14*(self.num_atoms() as f64)
        };
        self.error = new_error + self.error;
        if self.error > self.expected_accuracy(new_e) {
            self.error *= 0.0;
            self.E = self.compute_energy();
        } else {
            self.E = new_e;
        }
    }
    fn expected_accuracy(&self, newe: Energy) -> Energy {
        newe.abs()*1e-14*(self.num_atoms() as f64)*(self.num_atoms() as f64)
    }
}

impl From<WcaParams> for Wca {
    fn from(params: WcaParams) -> Wca {
        let cell = Cell::new(&params._dim, r_cutoff());
        if cell.r_cutoff > cell.box_diagonal.x ||
           cell.r_cutoff > cell.box_diagonal.y ||
           cell.r_cutoff > cell.box_diagonal.z
        {
            panic!("The cell is not large enough for the well width, sorry!");
        }
        Wca {
            E: 0.0*units::EPSILON,
            error: 0.0*units::EPSILON,
            cell,
            possible_change: Change::None,
        }
    }
}

/// This defines the pressure.
#[derive(Serialize, Deserialize, Debug, Default)]
pub struct Collected {
    /// The sum of dE/dϵ where ϵ is the strain (or fractional
    /// expansion) over all elements in the ensemble.
    pub pexc_tot: f64,
    /// The number times we have summed up the excess_pressure_total
    /// above.  Divide by this (and by the volume) to get the excess
    /// pressure.
    pub count: u64,
}

impl System for Wca {
    type CollectedData = Collected;
    fn collect_data(&self, data: &mut Collected, iter: u64) {
        if iter % ((self.num_atoms()*self.num_atoms()) as u64) == 0 {
            data.count += 1;
            let mut p: Energy = units::EPSILON*0.0;
            for (which, &r1) in self.cell.positions.iter().enumerate() {
                for r2 in self.cell.maybe_interacting_atoms_excluding(r1, which) {
                    p += potential_pressure((r1-r2).norm2());
                }
            }
            data.pexc_tot += *(p/units::EPSILON).value();
        }
    }
    fn energy(&self) -> Energy {
        self.E
    }
    fn compute_energy(&self) -> Energy {
        let mut e: Energy = units::EPSILON*0.0;
        for (which, &r1) in self.cell.positions.iter().enumerate() {
            for r2 in self.cell.maybe_interacting_atoms_excluding(r1, which) {
                e += potential((r1-r2).norm2());
            }
        }
        e*0.5
    }
    fn update_caches(&mut self) {
        self.cell.update_caches();
    }
    fn lowest_possible_energy(&self) -> Option<Energy> {
        Some(0.0*units::EPSILON)
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

impl ConfirmSystem for Wca {
    fn confirm(&mut self) {
        match self.possible_change {
            Change::None => (),
            Change::Move{which, to, e} => {
                self.cell.move_atom(which, to);
                self.set_energy(e);
                self.possible_change = Change::None;
            },
            Change::Add{to, e} => {
                self.cell.add_atom_at(to);
                self.set_energy(e);
                self.possible_change = Change::None;
            },
            Change::Remove{which, e} => {
                self.cell.remove_atom(which);
                self.set_energy(e);
                self.possible_change = Change::None;
            },
        }
    }
}

impl GrandSystem for Wca {
    fn plan_add(&mut self, rng: &mut MyRng) -> Option<Energy> {
        let r = self.cell.put_in_cell(
            Vector3d::new(Length::new(rng.sample(Uniform::new(0.0, self.cell.box_diagonal.x.value_unsafe))),
                          Length::new(rng.sample(Uniform::new(0.0, self.cell.box_diagonal.y.value_unsafe))),
                          Length::new(rng.sample(Uniform::new(0.0, self.cell.box_diagonal.z.value_unsafe)))));
        self.add_atom_at(r)
    }
    fn plan_remove(&mut self, rng: &mut MyRng) -> Energy {
        let which = rng.sample(Uniform::new(0, self.cell.positions.len()));
        self.remove_atom_number(which)
    }
    fn num_atoms(&self) -> usize {
        self.cell.positions.len()
    }
}

impl MovableSystem for Wca {
    fn plan_move(&mut self, rng: &mut MyRng, mean_distance: Length) -> Option<Energy> {
        use crate::rng::vector;
        if self.cell.positions.len() > 0 {
            let which = rng.sample(Uniform::new(0, self.cell.positions.len()));
            let to = self.cell.put_in_cell(unsafe { *self.cell.positions.get_unchecked(which) } + vector(rng)*mean_distance);
            self.move_atom(which, to)
        } else {
            None
        }
    }
    fn max_size(&self) -> Length {
        self.cell.box_diagonal.norm2().sqrt()
    }
}



/// A description of the cell dimensions and number.
#[derive(Serialize, Deserialize, Debug, ClapMe)]
pub enum CellDimensionsGivenNumber {
    /// The three widths of the cell
    CellWidth(Vector3d<Length>),
    /// The volume of the cell
    CellVolume(Volume),
    /// The reduced density, from which we can compute volume
    ReducedDensity(units::Density<f64>),
}

/// Parameters needed to configure a finite-N WCAl system.
#[derive(Serialize, Deserialize, Debug, ClapMe)]
#[allow(non_snake_case)]
pub struct WcaNParams {
    /// The sice of the cell.
    pub _dim: CellDimensionsGivenNumber,
    /// The number of atoms.
    pub N: usize,
}

impl Default for WcaNParams {
    fn default() -> Self {
        WcaNParams {
            _dim: CellDimensionsGivenNumber::ReducedDensity(units::Density::new(1.0)),
            N: 100,
        }
    }
}

impl From<WcaNParams> for Wca {
    fn from(params: WcaNParams) -> Wca {
        let n = params.N;
        let dim: CellDimensions = match params._dim {
            CellDimensionsGivenNumber::CellWidth(v)
                => CellDimensions::CellWidth(v),
            CellDimensionsGivenNumber::CellVolume(v)
                => CellDimensions::CellVolume(v),
            CellDimensionsGivenNumber::ReducedDensity(d)
                => CellDimensions::CellVolume((n as f64)/d),
        };
        let mut sw = Wca::from(WcaParams {
            _dim: dim,
        });

        // Atoms will be initially placed on a face centered cubic (fcc) grid
        // Note that the unit cells need not be actually "cubic", but the fcc grid will
        //   be stretched to cell dimensions
        let min_cell_width = 2.0*2.0f64.sqrt()*units::R; // minimum cell width
        let spots_per_cell = 4; // spots in each fcc periodic unit cell

        // cells holds the max number of cells that will fit in the x,
        // y, and z dimensions
        let cells = [*(sw.cell.box_diagonal.x/min_cell_width).value() as usize,
                     *(sw.cell.box_diagonal.y/min_cell_width).value() as usize,
                     *(sw.cell.box_diagonal.z/min_cell_width).value() as usize];
        // It is usefull to know our cell dimensions
        let cell_width = [sw.cell.box_diagonal.x/cells[0] as f64,
                          sw.cell.box_diagonal.y/cells[1] as f64,
                          sw.cell.box_diagonal.z/cells[2] as f64];
        for i in 0..3 {
            assert!(cell_width[i] >= min_cell_width);
        }
        // Define ball positions relative to cell position
        let offset = [Vector3d::new(Length::new(0.0), Length::new(0.0), Length::new(0.0)),
                      Vector3d::new(Length::new(0.0), cell_width[1], cell_width[2])/2.0,
                      Vector3d::new(cell_width[0], Length::new(0.0), cell_width[2])/2.0,
                      Vector3d::new(cell_width[0], cell_width[1], Length::new(0.0))/2.0];
        // Reserve total_spots at random to be occupied
        let total_spots = spots_per_cell*cells[0]*cells[1]*cells[1];
        if total_spots < params.N {
            println!("problem:  {} < {}", total_spots, params.N);
        }
        assert!(total_spots >= params.N);
        let mut spots_reserved = vec![vec![vec![[false; 4]; cells[2]]; cells[1]]; cells[0]];
        let mut rng = crate::rng::MyRng::seed_from_u64(0);
        for _ in 0..params.N {
            loop {
                // This is an inefficient but relatively
                // hard-to-get-wrong way to randomly sample spots.
                // Speed shouldn't matter here.
                let i = rng.sample(Uniform::new(0, cells[0]));
                let j = rng.sample(Uniform::new(0, cells[1]));
                let k = rng.sample(Uniform::new(0, cells[2]));
                let l = rng.sample(Uniform::new(0, 4));
                if !spots_reserved[i][j][k][l] {
                    spots_reserved[i][j][k][l] = true;
                    sw.add_atom_at(Vector3d::new(i as f64*cell_width[0],
                                                 j as f64*cell_width[1],
                                                 k as f64*cell_width[2])
                                   + offset[l]);
                    sw.confirm();
                    sw.verify_energy();
                    // assert_eq!(sw.energy(), sw.compute_energy());
                    break;
                }
            }
        }
        sw
    }
}

#[cfg(test)]
fn closest_distance_matches(natoms: usize) {
    let mut sw = mk_wca(natoms, 0.3);
    for &r1 in sw.cell.positions.iter() {
        for &r2 in sw.cell.positions.iter() {
            assert_eq!(sw.cell.closest_distance2(r1,r2),
                       sw.cell.unsafe_closest_distance2(r1,r2));
            assert_eq!(sw.cell.closest_distance2(r1,r2),
                       sw.cell.sloppy_closest_distance2(r1,r2));
        }
    }
    let mut rng = MyRng::seed_from_u64(1);
    for _ in 0..100000 {
        sw.plan_move(&mut rng, Length::new(1.0));
        sw.confirm();
    }
    for &r1 in sw.cell.positions.iter() {
        for &r2 in sw.cell.positions.iter() {
            assert_eq!(sw.cell.closest_distance2(r1,r2),
                       sw.cell.unsafe_closest_distance2(r1,r2));
            assert_eq!(sw.cell.closest_distance2(r1,r2),
                       sw.cell.sloppy_closest_distance2(r1,r2));
        }
    }
}

#[test]
fn closest_distance_matches_n50() {
    closest_distance_matches(50);
}
#[test]
fn closest_distance_matches_n100() {
    closest_distance_matches(100);
}
#[test]
fn closest_distance_matches_n200() {
    closest_distance_matches(200);
}

#[cfg(test)]
fn maybe_interacting_needs_no_shifting(natoms: usize) {
    let mut sw = mk_wca(natoms, 0.3);
    let mut rng = MyRng::seed_from_u64(1);
    for _ in 0..100000 {
        sw.plan_move(&mut rng, Length::new(1.0));
        sw.confirm();
    }
    for &r1 in sw.cell.positions.iter() {
        for r2 in sw.cell.maybe_interacting_atoms(r1) {
            println!("comparing positions:\n{}\n  and\n{}", r1, r2);
            assert_eq!(sw.cell.closest_distance2(r1,r2),
                       sw.cell.unsafe_closest_distance2(r1,r2));
            assert_eq!(sw.cell.closest_distance2(r1,r2),
                       sw.cell.sloppy_closest_distance2(r1,r2));
            assert_eq!(sw.cell.closest_distance2(r1,r2),
                       (r1-r2).norm2());
        }
    }
    println!("starting with exclusion.");
    for (which, &r1) in sw.cell.positions.iter().enumerate() {
        for r2 in sw.cell.maybe_interacting_atoms_excluding(r1, which) {
            assert_eq!(sw.cell.closest_distance2(r1,r2),
                       sw.cell.unsafe_closest_distance2(r1,r2));
            assert_eq!(sw.cell.closest_distance2(r1,r2),
                       sw.cell.sloppy_closest_distance2(r1,r2));
            assert_eq!(sw.cell.closest_distance2(r1,r2),
                       (r1-r2).norm2());
        }
    }
}

// #[test]
// fn maybe_interacting_needs_no_shifting_n50() {
//     maybe_interacting_needs_no_shifting(50);
// }
#[test]
fn maybe_interacting_needs_no_shifting_n100() {
    maybe_interacting_needs_no_shifting(100);
}
#[test]
fn maybe_interacting_needs_no_shifting_n200() {
    maybe_interacting_needs_no_shifting(200);
}

#[test]
fn maybe_interacting_includes_everything() {
    let mut sw = mk_wca(100, 0.3);
    let mut rng = MyRng::seed_from_u64(1);
    for _ in 0..100000 {
        sw.plan_move(&mut rng, Length::new(1.0));
        sw.confirm()
    }
    for &r1 in sw.cell.positions.iter() {
        sw.cell.verify_maybe_interacting_includes_everything(r1);
    }
}

#[cfg(test)]
fn maybe_interacting_excluding_includes_everything(natoms: usize) {
    let mut sw = mk_wca(natoms, 0.3);
    let mut rng = MyRng::seed_from_u64(1);
    for (which, &r1) in sw.cell.positions.iter().enumerate() {
        sw.cell.verify_maybe_interacting_excluding_includes_everything(r1, which);
    }
    for _ in 0..100000 {
        sw.plan_move(&mut rng, Length::new(1.0));
        sw.confirm();
    }
    println!("Finished moving stuff around...");
    for (which, &r1) in sw.cell.positions.iter().enumerate() {
        sw.cell.verify_maybe_interacting_excluding_includes_everything(r1, which);
        sw.cell.verify_maybe_interacting_excluding_includes_everything(r1, 0);
    }
}

#[test]
fn maybe_interacting_excluding_includes_everything_n50() {
    maybe_interacting_excluding_includes_everything(50);
}

#[test]
fn maybe_interacting_excluding_includes_everything_n100() {
    maybe_interacting_excluding_includes_everything(100);
}

#[test]
fn maybe_interacting_excluding_includes_everything_n200() {
    maybe_interacting_excluding_includes_everything(200);
}

#[test]
fn energy_is_right_n3() {
    energy_is_right(3, 0.1);
}
#[test]
fn energy_is_right_n50() {
    energy_is_right(50, 0.3);
}
#[test]
fn energy_is_right_n100() {
    energy_is_right(100, 0.3);
}
#[test]
fn energy_is_right_n200() {
    energy_is_right(200, 0.3);
}

#[cfg(test)]
fn energy_is_right(natoms: usize, ff: f64) {
    let mut sw = mk_wca(natoms, ff);
    sw.verify_energy();
    let mut rng = MyRng::seed_from_u64(1);
    let mut old_energy = sw.energy();
    let maxe = (natoms as f64)*16.0*units::EPSILON;
    let mut i = 0.0;
    while i < 1000.0 {
        if let Some(newe) = sw.plan_move(&mut rng, Length::new(1.0)) {
            if newe < maxe || newe < old_energy {
                sw.confirm();
                println!("after move {}... {} vs {}", i, sw.energy(), sw.compute_energy());
                sw.verify_energy();
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
fn mk_wca(natoms: usize, n: f64) -> Wca {
    let mut param = WcaNParams::default();
    param._dim = CellDimensionsGivenNumber::ReducedDensity(units::Density::new(n));
    param.N = natoms;
    Wca::from(param)
}
