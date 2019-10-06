//! A square well fluid.

use super::*;

use dimensioned::Dimensionless;
use dimensioned::{Cbrt,Abs};
use vector3d::Vector3d;
use rand::prelude::*;
use rand::distributions::Uniform;
use std::f64::consts::PI;
use std::default::Default;

/// A description of the cell dimensions
#[derive(Serialize, Deserialize, Debug, ClapMe)]
#[allow(non_snake_case)]
pub enum CellDimensions {
    /// The three widths of the cell
    CellWidth(Vector3d<Length>),
    /// The volume of the cell
    CellVolume(Volume),
}

/// The parameters needed to configure a square well system.
#[derive(Serialize, Deserialize, Debug, ClapMe)]
pub struct SquareWellParams {
    well_width: Unitless,
    _dim: CellDimensions,
}

#[allow(non_snake_case)]
/// A square well fluid.
#[derive(Serialize, Deserialize, Debug)]
pub struct SquareWell {
    /// The dimensionless well width.
    well_width: Length,
    /// The atom positions
    pub positions: Vec<Vector3d<Length>>,
    /// The energy of the system
    E: Energy,
    /// The dimensions of the box.
    box_diagonal: Vector3d<Length>,
    /// The last change we made (and might want to undo).
    possible_change: Change,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
enum Change {
    Move { which: usize, to: Vector3d<Length>, e: Energy },
    /// The added atom is always pushed to the end of the vector!
    Add { to: Vector3d<Length>, e: Energy },
    Remove { which: usize, e: Energy },
    None,
}

impl SquareWell {
    fn max_interaction(&self) -> u64 {
        max_balls_within(self.well_width)
    }
    /// Add an atom at a given location.  Returns the change in
    /// energy, or `None` if the atom could not be placed there.
    pub fn add_atom_at(&mut self, r: Vector3d<Length>) -> Option<Energy> {
        let mut e = self.E;
        for &r1 in self.positions.iter() {
            let dist2 = self.closest_distance2(r1,r);
            if dist2 < units::SIGMA*units::SIGMA {
                self.possible_change = Change::None;
                return None;
            } else if dist2 < self.well_width*self.well_width {
                e -= units::EPSILON;
            }
        }
        self.possible_change = Change::Add{ to: r, e };
        Some(e)
    }
    /// Plan to move a specified atom.  Returns the energy, or `None`
    /// if the atom could not be placed there.
    pub fn move_atom(&mut self, which: usize, r: Vector3d<Length>) -> Option<Energy> {
        let mut e = self.E;
        let from = self.positions[which];
        for &r1 in self.positions.iter() {
            let old_dist2 = self.closest_distance2(r1,from);
            if old_dist2 > Area::new(0.0) { // not the same atom!
                let dist2 = self.closest_distance2(r1,r);
                if dist2 < units::SIGMA*units::SIGMA {
                    self.possible_change = Change::None;
                    return None;
                } else if dist2 < self.well_width*self.well_width {
                    e -= units::EPSILON;
                }
                if old_dist2 < self.well_width*self.well_width {
                    e += units::EPSILON;
                }
            }
        }
        self.possible_change = Change::Move{ which, to: r, e };
        Some(e)
    }
    /// Plan to remove the specified atom.  Returns the energy.
    pub fn remove_atom_number(&mut self, which: usize) -> Energy {
        let r = self.positions[which];
        let mut e = self.E;
        for &r1 in self.positions.iter().filter(|&&r1| r1 != r) {
            if self.closest_distance2(r1,r) < self.well_width*self.well_width {
                e += units::EPSILON;
            }
        }
        self.possible_change = Change::Remove{ which, e };
        e
    }
    /// PUBLIC FOR TESTING ONLY! The shortest distance squared between two vectors.
    pub fn closest_distance2(&self, r1: Vector3d<Length>, r2: Vector3d<Length>) -> Area {
        let mut dr = r2 - r1;
        if dr.x < -0.5*self.box_diagonal.x {
            while {
                dr.x += self.box_diagonal.x;
                dr.x < -0.5*self.box_diagonal.x
            } {}
        } else {
            while dr.x > 0.5*self.box_diagonal.x {
                dr.x -= self.box_diagonal.x;
            }
        }
        if dr.y < -0.5*self.box_diagonal.y {
            while {
                dr.y += self.box_diagonal.y;
                dr.y < -0.5*self.box_diagonal.y
            } {}
        } else {
            while dr.y > 0.5*self.box_diagonal.y {
                dr.y -= self.box_diagonal.y;
            }
        }
        if dr.z < -0.5*self.box_diagonal.z {
            while {
                dr.z += self.box_diagonal.z;
                dr.z < -0.5*self.box_diagonal.z
            } {}
        } else {
            while dr.z > 0.5*self.box_diagonal.z {
                dr.z -= self.box_diagonal.z;
            }
        }
        dr.norm2()
    }
    /// PUBLIC FOR TESTING ONLY! The shortest distance squared between two vectors.
    pub fn sloppy_closest_distance2(&self, r1: Vector3d<Length>, r2: Vector3d<Length>) -> Area {
        let mut dr = r2 - r1;
        while dr.x < -0.5*self.box_diagonal.x {
            dr.x += self.box_diagonal.x;
        }
        while dr.x > 0.5*self.box_diagonal.x {
            dr.x -= self.box_diagonal.x;
        }
        while dr.y < -0.5*self.box_diagonal.y {
            dr.y += self.box_diagonal.y;
        }
        while dr.y > 0.5*self.box_diagonal.y {
            dr.y -= self.box_diagonal.y;
        }
        while dr.z < -0.5*self.box_diagonal.z {
            dr.z += self.box_diagonal.z;
        }
        while dr.z > 0.5*self.box_diagonal.z {
            dr.z -= self.box_diagonal.z;
        }
        dr.norm2()
    }
    /// PUBLIC FOR TESTING ONLY! The shortest distance squared between two vectors.
    pub fn unsafe_closest_distance2(&self, r1: Vector3d<Length>, r2: Vector3d<Length>) -> Area {
        let mut dr = r2 - r1;
        if dr.x < -0.5*self.box_diagonal.x {
            dr.x += self.box_diagonal.x;
        } else if dr.x > 0.5*self.box_diagonal.x {
            dr.x -= self.box_diagonal.x;
        }
        if dr.y < -0.5*self.box_diagonal.y {
            dr.y += self.box_diagonal.y;
        } else if dr.y > 0.5*self.box_diagonal.y {
            dr.y -= self.box_diagonal.y;
        }
        if dr.z < -0.5*self.box_diagonal.z {
            dr.z += self.box_diagonal.z;
        } else if dr.z > 0.5*self.box_diagonal.z {
            dr.z -= self.box_diagonal.z;
        }
        dr.norm2()
    }
    /// PUBLIC FOR TESTING ONLY! Put position into the cell.
    pub fn put_in_cell(&self, mut r: Vector3d<Length>) -> Vector3d<Length> {
        if r.x < Length::new(0.0) {
            while {
                r.x += self.box_diagonal.x;
                r.x < Length::new(0.0)
            } {}
        } else {
            while r.x > self.box_diagonal.x {
                r.x -= self.box_diagonal.x;
            }
        }
        if r.y < Length::new(0.0) {
            while {
                r.y += self.box_diagonal.y;
                r.y < Length::new(0.0)
            } {}
        } else {
            while r.y > self.box_diagonal.y {
                r.y -= self.box_diagonal.y;
            }
        }
        if r.z < Length::new(0.0) {
            while {
                r.z += self.box_diagonal.z;
                r.z < Length::new(0.0)
            } {}
        } else {
            while r.z > self.box_diagonal.z {
                r.z -= self.box_diagonal.z;
            }
        }
        r
    }
    /// PUBLIC FOR TESTING ONLY! Put position into the cell.
    pub fn sloppy_put_in_cell(&self, mut r: Vector3d<Length>) -> Vector3d<Length> {
        while r.x < Length::new(0.0) {
            r.x += self.box_diagonal.x;
        }
        while r.x > self.box_diagonal.x {
            r.x -= self.box_diagonal.x;
        }
        while r.y < Length::new(0.0) {
            r.y += self.box_diagonal.y;
        }
        while r.y > self.box_diagonal.y {
            r.y -= self.box_diagonal.y;
        }
        while r.z < Length::new(0.0) {
            r.z += self.box_diagonal.z;
        }
        while r.z > self.box_diagonal.z {
            r.z -= self.box_diagonal.z;
        }
        r
    }
}

impl From<SquareWellParams> for SquareWell {
    fn from(params: SquareWellParams) -> SquareWell {
        let box_diagonal = match params._dim {
            CellDimensions::CellWidth(w) => {
                Vector3d::new(w.x.abs(),w.y.abs(),w.z.abs())
            },
            CellDimensions::CellVolume(v) => {
                let w = v.cbrt();
                Vector3d::new(w,w,w)
            }
        };
        SquareWell {
            well_width: params.well_width*units::SIGMA,
            positions: Vec::new(),
            E: 0.0*units::EPSILON,
            box_diagonal: box_diagonal,
            possible_change: Change::None,
        }
    }
}

impl System for SquareWell {
    type CollectedData = ();
    fn energy(&self) -> Energy {
        self.E
    }
    fn compute_energy(&self) -> Energy {
        let mut e: Energy = units::EPSILON*0.0;
        for (i, &r1) in self.positions[..self.positions.len()-1].iter().enumerate() {
            for &r2 in self.positions[i+1..].iter() {
                if self.closest_distance2(r1,r2) < self.well_width*self.well_width {
                    e -= units::EPSILON;
                }
            }
        }
        e
    }
    fn delta_energy(&self) -> Option<Energy> {
        Some(units::EPSILON)
    }
    fn greatest_possible_energy(&self) -> Option<Energy> {
        Some(0.0*units::EPSILON)
    }
    fn lowest_possible_energy(&self) -> Option<Energy> {
        Some(-(self.positions.len() as f64)*(self.max_interaction() as f64)*units::EPSILON)
    }
}

impl ConfirmSystem for SquareWell {
    fn confirm(&mut self) {
        match self.possible_change {
            Change::None => (),
            Change::Move{which, to, e} => {
                self.positions[which] = to;
                self.E = e;
                self.possible_change = Change::None;
            },
            Change::Add{to, e} => {
                self.positions.push(to);
                self.E = e;
                self.possible_change = Change::None;
            },
            Change::Remove{which, e} => {
                self.positions.swap_remove(which);
                self.E = e;
                self.possible_change = Change::None;
            },
        }
    }
}

impl GrandSystem for SquareWell {
    fn plan_add(&mut self, rng: &mut MyRng) -> Option<Energy> {
        let r = self.put_in_cell(Vector3d::new(Length::new(rng.sample(Uniform::new(0.0, self.box_diagonal.x.value_unsafe))),
                                               Length::new(rng.sample(Uniform::new(0.0, self.box_diagonal.y.value_unsafe))),
                                               Length::new(rng.sample(Uniform::new(0.0, self.box_diagonal.z.value_unsafe)))));
        self.add_atom_at(r)
    }
    fn plan_remove(&mut self, rng: &mut MyRng) -> Energy {
        let which = rng.sample(Uniform::new(0, self.positions.len()));
        self.remove_atom_number(which)
    }
    fn num_atoms(&self) -> usize {
        self.positions.len()
    }
}

impl MovableSystem for SquareWell {
    fn plan_move(&mut self, rng: &mut MyRng, mean_distance: Length) -> Option<Energy> {
        if self.positions.len() > 0 {
            let which = rng.sample(Uniform::new(0, self.positions.len()));
            let to = self.put_in_cell(self.positions[which] + rng.vector()*mean_distance);
            self.move_atom(which, to)
        } else {
            None
        }
    }
}

fn max_balls_within(mut distance: Length) -> u64 {
    distance += 1e-10*units::SIGMA; // add a tiny, but necessary margin of error
    let a = 2.0_f64.sqrt()*units::SIGMA; // fcc lattice constant
    let c = (*(distance/a).value()).ceil() as i64 + 1; // number of cubic fcc cells to go out from center
    let mut num: i64 = -1; // number of balls within a given radius; don't count the center ball
    let d2 = distance*distance;

    for n in -c .. c+1 {
        for m in -c .. c+1 {
            for l in -c .. c+1 {
                let x0 = ((m+l) as f64)*a;
                let y0 = ((n+l) as f64)*a;
                let z0 = ((m+n) as f64)*a;
                if x0*x0 + y0*y0 + z0*z0 <= d2 {
                    num += 1;
                }
                if (x0 + 0.5*a)*(x0 + 0.5*a) + (y0 + 0.5*a)*(y0 + 0.5*a) + z0*z0 <= d2 {
                    num += 1;
                }
                if (x0 + 0.5*a)*(x0 + 0.5*a) + y0*y0 + (z0 + 0.5*a)*(z0 + 0.5*a) <= d2 {
                    num += 1;
                }
                if x0*x0 + (y0 + 0.5*a)*(y0 + 0.5*a) + (z0 + 0.5*a)*(z0 + 0.5*a) <= d2 {
                    num += 1;
                }
            }
        }
    }
    num as u64
}



/// A description of the cell dimensions and number.
#[derive(Serialize, Deserialize, Debug, ClapMe)]
pub enum CellDimensionsGivenNumber {
    /// The three widths of the cell
    CellWidth(Vector3d<Length>),
    /// The volume of the cell
    CellVolume(Volume),
    /// The filling fraction, from which we compute volume
    FillingFraction(Unitless),
}

/// Parameters needed to configure a finite-N square-well system.
#[derive(Serialize, Deserialize, Debug, ClapMe)]
#[allow(non_snake_case)]
pub struct SquareWellNParams {
    /// The width of the well, relative to the diameter.
    pub well_width: Unitless,
    /// The sice of the cell.
    pub _dim: CellDimensionsGivenNumber,
    /// The number of atoms.
    pub N: usize,
}

impl Default for SquareWellNParams {
    fn default() -> Self {
        SquareWellNParams {
            well_width: Unitless::new(1.3),
            _dim: CellDimensionsGivenNumber::FillingFraction(Unitless::new(0.3)),
            N: 100,
        }
    }
}

impl From<SquareWellNParams> for SquareWell {
    fn from(params: SquareWellNParams) -> SquareWell {
        let n = params.N;
        let dim: CellDimensions = match params._dim {
            CellDimensionsGivenNumber::CellWidth(v)
                => CellDimensions::CellWidth(v),
            CellDimensionsGivenNumber::CellVolume(v)
                => CellDimensions::CellVolume(v),
            CellDimensionsGivenNumber::FillingFraction(f)
                => CellDimensions::CellVolume((n as f64)*(PI*units::SIGMA*units::SIGMA*units::SIGMA/6.0)/f),
        };
        let mut sw = SquareWell::from(SquareWellParams {
            _dim: dim,
            well_width: params.well_width,
        });

        // Atoms will be initially placed on a face centered cubic (fcc) grid
        // Note that the unit cells need not be actually "cubic", but the fcc grid will
        //   be stretched to cell dimensions
        let min_cell_width = 2.0*2.0f64.sqrt()*units::R; // minimum cell width
        let spots_per_cell = 4; // spots in each fcc periodic unit cell

        // cells holds the max number of cells that will fit in the x,
        // y, and z dimensions
        let cells = [*(sw.box_diagonal.x/min_cell_width).value() as usize,
                     *(sw.box_diagonal.y/min_cell_width).value() as usize,
                     *(sw.box_diagonal.z/min_cell_width).value() as usize];
        // It is usefull to know our cell dimensions
        let cell_width = [sw.box_diagonal.x/cells[0] as f64,
                          sw.box_diagonal.y/cells[1] as f64,
                          sw.box_diagonal.z/cells[2] as f64];
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
        let mut rng = ::rng::MyRng::from_u64(0);
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
                    break;
                }
            }
        }
        sw
    }
}

#[test]
fn closest_distance_matches() {
    use std::default::Default;
    let mut sw = SquareWell::from(SquareWellNParams::default());
    for &r1 in sw.positions.iter() {
        for &r2 in sw.positions.iter() {
            assert_eq!(sw.closest_distance2(r1,r2), sw.unsafe_closest_distance2(r1,r2));
            assert_eq!(sw.closest_distance2(r1,r2), sw.sloppy_closest_distance2(r1,r2));
        }
    }
    let mut rng = MyRng::from_u64(1);
    for _ in 0..100000 {
        sw.plan_move(&mut rng, Length::new(1.0));
        sw.confirm();
    }
    for &r1 in sw.positions.iter() {
        for &r2 in sw.positions.iter() {
            assert_eq!(sw.closest_distance2(r1,r2), sw.unsafe_closest_distance2(r1,r2));
            assert_eq!(sw.closest_distance2(r1,r2), sw.sloppy_closest_distance2(r1,r2));
        }
    }
}

#[test]
fn energy_is_right() {
    use std::default::Default;
    let mut sw = SquareWell::from(SquareWellNParams::default());
    assert_eq!(sw.energy(), sw.compute_energy());
    let mut rng = MyRng::from_u64(1);
    for _ in 0..100 {
        println!("making a move...");
        sw.plan_move(&mut rng, Length::new(1.0));
        sw.confirm();
        assert_eq!(sw.energy(), sw.compute_energy());
    }
}
