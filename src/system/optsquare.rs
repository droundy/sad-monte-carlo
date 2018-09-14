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

/// A square well fluid.
#[derive(Serialize, Deserialize, Debug)]
pub struct Cell {
    /// The dimensions of the box.
    box_diagonal: Vector3d<Length>,
    /// The dimensionless well width.
    well_width: Length,
    /// The atom positions
    pub positions: Vec<Vector3d<Length>>,
    /// The dimensions of the subcell grid
    num_subcells: Vector3d<usize>,
    /// The subcell lists
    subcells: Vec<Vec<Vector3d<Length>>>
}

#[allow(non_snake_case)]
/// A square well fluid.
#[derive(Serialize, Deserialize, Debug)]
pub struct SquareWell {
    /// The energy of the system
    E: Energy,
    /// The dimensions of the box.
    pub cell: Cell,
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

fn modulus(i: isize, sz: usize) -> usize {
    let sz = sz as isize;
    (((i % sz) + sz) % sz) as usize
}

impl ::std::ops::Index<Vector3d<isize>> for Cell {
    type Output = Vec<Vector3d<Length>>;
    fn index(&self, i: Vector3d<isize>) -> &Vec<Vector3d<Length>> {
        let ix = modulus(i.x, self.num_subcells.x);
        let iy = modulus(i.y, self.num_subcells.y);
        let iz = modulus(i.z, self.num_subcells.z);
        // &self.subcells[ix*(self.num_subcells.y*self.num_subcells.z) + iy*self.num_subcells.z + iz]
        unsafe { self.subcells.get_unchecked(ix*(self.num_subcells.y*self.num_subcells.z) + iy*self.num_subcells.z + iz) }
    }
}
impl ::std::ops::IndexMut<Vector3d<isize>> for Cell {
    fn index_mut(&mut self, i: Vector3d<isize>) -> &mut Vec<Vector3d<Length>> {
        let ix = modulus(i.x, self.num_subcells.x);
        let iy = modulus(i.y, self.num_subcells.y);
        let iz = modulus(i.z, self.num_subcells.z);
        // &mut self.subcells[ix*(self.num_subcells.y*self.num_subcells.z) + iy*self.num_subcells.z + iz]
        unsafe { self.subcells.get_unchecked_mut(ix*(self.num_subcells.y*self.num_subcells.z) + iy*self.num_subcells.z + iz) }
    }
}

const NO_NEIGHBORS: [Vector3d<isize>; 0] = [];
const NEIGHBORS: [Vector3d<isize>; 26] = [
    Vector3d { x:  1, y:  0, z:  0 },
    Vector3d { x: -1, y:  0, z:  0 },
    Vector3d { x:  0, y:  1, z:  0 },
    Vector3d { x:  0, y: -1, z:  0 },
    Vector3d { x:  0, y:  0, z:  1 },
    Vector3d { x:  0, y:  0, z: -1 },
    Vector3d { x:  0, y:  1, z:  1 },
    Vector3d { x:  0, y:  1, z: -1 },
    Vector3d { x:  0, y: -1, z:  1 },
    Vector3d { x:  0, y: -1, z: -1 },
    Vector3d { x:  1, y:  0, z:  1 },
    Vector3d { x:  1, y:  0, z: -1 },
    Vector3d { x: -1, y:  0, z:  1 },
    Vector3d { x: -1, y:  0, z: -1 },
    Vector3d { x:  1, y:  1, z:  0 },
    Vector3d { x:  1, y: -1, z:  0 },
    Vector3d { x: -1, y:  1, z:  0 },
    Vector3d { x: -1, y: -1, z:  0 },
    Vector3d { x:  1, y:  1, z:  1 },
    Vector3d { x: -1, y:  1, z:  1 },
    Vector3d { x:  1, y: -1, z:  1 },
    Vector3d { x:  1, y:  1, z: -1 },
    Vector3d { x:  1, y: -1, z: -1 },
    Vector3d { x: -1, y:  1, z: -1 },
    Vector3d { x: -1, y: -1, z:  1 },
    Vector3d { x: -1, y: -1, z: -1 }];

struct NewNeighborIterator<'a> {
    cell: &'a Cell,
    sc: Vector3d<isize>,
    offset_iter: ::std::slice::Iter<'a, Vector3d<isize>>,
    position_iter: ::std::slice::Iter<'a, Vector3d<Length>>,
    shift: Vector3d<Length>,
    exclude: Option<Vector3d<Length>>,
    shift_near: Option<Vector3d<Length>>,
}
impl<'a> Iterator for NewNeighborIterator<'a> {
    type Item = Vector3d<Length>;
    fn next(&mut self) -> Option<Vector3d<Length>> {
        loop {
            while let Some(&r) = self.position_iter.next() {
                if Some(r) != self.exclude {
                    if let Some(shift_near) = self.shift_near {
                        let mut me = r;
                        if me.z - shift_near.z < -self.cell.box_diagonal.z*0.5 {
                            me.z += self.cell.box_diagonal.z;
                        } else if me.z - shift_near.z > self.cell.box_diagonal.z*0.5 {
                            me.z -= self.cell.box_diagonal.z;
                        }
                        if me.y - shift_near.y < -self.cell.box_diagonal.y*0.5 {
                            me.y += self.cell.box_diagonal.y;
                        } else if me.y - shift_near.y > self.cell.box_diagonal.y*0.5 {
                            me.y -= self.cell.box_diagonal.y;
                        }
                        if me.x - shift_near.x < -self.cell.box_diagonal.x*0.5 {
                            me.x += self.cell.box_diagonal.x;
                        } else if me.x - shift_near.x > self.cell.box_diagonal.x*0.5 {
                            me.x -= self.cell.box_diagonal.x;
                        }
                        return Some(me);
                    } else {
                        return Some(r + self.shift);
                    }
                }
            }
            if let Some(&offset) = self.offset_iter.next() {
                let pos = &self.cell[self.sc + offset];
                self.position_iter = pos.iter();
                if pos.len() > 0 {
                self.shift = Vector3d::new(
                    if self.sc.x + offset.x < 0 {
                        -self.cell.box_diagonal.x
                    } else if self.sc.x + offset.x >= self.cell.num_subcells.x as isize {
                        self.cell.box_diagonal.x
                    } else {
                        Length::new(0.)
                    },
                    if self.sc.y + offset.y < 0 {
                        -self.cell.box_diagonal.y
                    } else if self.sc.y + offset.y >= self.cell.num_subcells.y as isize {
                        self.cell.box_diagonal.y
                    } else {
                        Length::new(0.)
                    },
                    if self.sc.z + offset.z < 0 {
                        -self.cell.box_diagonal.z
                    } else if self.sc.z + offset.z >= self.cell.num_subcells.z as isize {
                        self.cell.box_diagonal.z
                    } else {
                        Length::new(0.)
                    });
                }
            } else {
                return None;
            }
        }
    }
}

impl Cell {
    /// Create a new Cell with specified dimensions and interaction length.
    pub fn new(dim: &CellDimensions, interaction_length: Length) -> Self {
        let box_diagonal = match dim {
            CellDimensions::CellWidth(w) => {
                Vector3d::new(w.x.abs(),w.y.abs(),w.z.abs())
            },
            CellDimensions::CellVolume(v) => {
                let w = v.cbrt();
                Vector3d::new(w,w,w)
            }
        };
        let mut cells_x = (box_diagonal.x/interaction_length).value().floor() as usize;
        let mut cells_y = (box_diagonal.y/interaction_length).value().floor() as usize;
        let mut cells_z = (box_diagonal.z/interaction_length).value().floor() as usize;
        if cells_z < 4 || cells_y < 4 || cells_x < 4 {
            cells_z = 1;
            cells_y = 1;
            cells_x = 1;
        }
        Cell {
            box_diagonal: box_diagonal,
            well_width: interaction_length,
            positions: Vec::new(),
            num_subcells: Vector3d::new(cells_x,cells_y,cells_z),
            subcells: vec![Vec::new(); cells_x*cells_y*cells_z],
        }
    }
    /// Atoms that may be within well_width of r.
    pub fn maybe_interacting_atoms<'a>(&'a self, r: Vector3d<Length>)
                                       -> impl Iterator<Item=Vector3d<Length>> + 'a {
        if self.num_subcells.x == 1 {
            NewNeighborIterator {
                cell: self,
                sc: Vector3d::new(0,0,0),
                offset_iter: NO_NEIGHBORS.iter(),
                position_iter: self.positions.iter(),
                shift: Vector3d::new(Length::new(0.),Length::new(0.),Length::new(0.)),
                exclude: None,
                shift_near: Some(r),
            }
        } else {
            let sc = self.get_subcell(r);
            NewNeighborIterator {
                cell: self,
                sc: sc,
                offset_iter: NEIGHBORS.iter(),
                position_iter: self[sc].iter(),
                shift: Vector3d::new(Length::new(0.),Length::new(0.),Length::new(0.)),
                exclude: None,
                shift_near: None,
            }
        }
    }
    /// Atoms that may be within well_width of r.  This excludes any
    /// atom that is located precisely at rprime.
    pub fn maybe_interacting_atoms_excluding<'a>(&'a self, r: Vector3d<Length>, rprime: Vector3d<Length>)
                                   -> impl Iterator<Item=Vector3d<Length>> + 'a {
        if self.num_subcells.x == 1 {
            NewNeighborIterator {
                cell: self,
                sc: Vector3d::new(0,0,0),
                offset_iter: NO_NEIGHBORS.iter(),
                position_iter: self.positions.iter(),
                shift: Vector3d::new(Length::new(0.),Length::new(0.),Length::new(0.)),
                exclude: Some(rprime),
                shift_near: Some(r),
            }
        } else {
            let sc = self.get_subcell(r);
            NewNeighborIterator {
                cell: self,
                sc: sc,
                offset_iter: NEIGHBORS.iter(),
                position_iter: self[sc].iter(),
                shift: Vector3d::new(Length::new(0.),Length::new(0.),Length::new(0.)),
                exclude: Some(rprime),
                shift_near: None,
            }
        }
    }
    /// Find the cell for a given vector.
    pub fn get_subcell(&self, r: Vector3d<Length>) -> Vector3d<isize> {
        Vector3d::new((r.x/self.box_diagonal.x*self.num_subcells.x as f64).value().floor() as isize,
                      (r.y/self.box_diagonal.y*self.num_subcells.y as f64).value().floor() as isize,
                      (r.z/self.box_diagonal.z*self.num_subcells.z as f64).value().floor() as isize)
    }
    /// Add an atom at a given location.
    pub fn add_atom_at(&mut self, r: Vector3d<Length>) {
        self.positions.push(r);
        let sc = self.get_subcell(r);
        self[sc].push(r);
    }
    /// Move the atom.  Return the previous position.
    pub fn move_atom(&mut self, which: usize, r: Vector3d<Length>) -> Vector3d<Length> {
        // let old = ::std::mem::replace(&mut self.positions[which], r);
        let old = ::std::mem::replace(unsafe { self.positions.get_unchecked_mut(which) }, r);
        let sc = self.get_subcell(r);
        self[sc].push(r);
        let oldsc = self.get_subcell(old);
        let mut oldi = 0;
        let subs = &mut self[oldsc];
        for (i,&v) in subs.iter().enumerate() {
            if v == old {
                oldi = i;
                break;
            }
        }
        subs.swap_remove(oldi);
        old
    }
    /// Remove an atom.  Return the previous position.
    pub fn remove_atom(&mut self, which: usize) -> Vector3d<Length> {
        let old = self.positions.swap_remove(which);
        let oldsc = self.get_subcell(old);
        let mut oldi = 0;
        for (i,&v) in self[oldsc].iter().enumerate() {
            if v == old {
                oldi = i;
                break;
            }
        }
        self[oldsc].swap_remove(oldi);
        old
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
            while r.x >= self.box_diagonal.x {
                r.x -= self.box_diagonal.x;
            }
        }
        if r.y < Length::new(0.0) {
            while {
                r.y += self.box_diagonal.y;
                r.y < Length::new(0.0)
            } {}
        } else {
            while r.y >= self.box_diagonal.y {
                r.y -= self.box_diagonal.y;
            }
        }
        if r.z < Length::new(0.0) {
            while {
                r.z += self.box_diagonal.z;
                r.z < Length::new(0.0)
            } {}
        } else {
            while r.z >= self.box_diagonal.z {
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
        while r.x >= self.box_diagonal.x {
            r.x -= self.box_diagonal.x;
        }
        while r.y < Length::new(0.0) {
            r.y += self.box_diagonal.y;
        }
        while r.y >= self.box_diagonal.y {
            r.y -= self.box_diagonal.y;
        }
        while r.z < Length::new(0.0) {
            r.z += self.box_diagonal.z;
        }
        while r.z >= self.box_diagonal.z {
            r.z -= self.box_diagonal.z;
        }
        r
    }
    /// for testing
    pub fn verify_subcells_include_everything(&self) {
        for &r in self.positions.iter() {
            let sc = self.get_subcell(r);
            assert!(self[sc].contains(&r));
        }
    }

    /// for testing
    pub fn verify_subcells_include_nothing_out_of_place(&self) {
        for ix in 0..self.num_subcells.x {
            for iy in 0..self.num_subcells.y {
                for iz in 0..self.num_subcells.z {
                    let sc = Vector3d::new(ix as isize,iy as isize,iz as isize);
                    for &r in self[sc].iter() {
                        assert_eq!(sc, self.get_subcell(r));
                    }
                }
            }
        }
        for &r in self.positions.iter() {
            let sc = self.get_subcell(r);
            if !self[sc].contains(&r) {
                println!("\n\nWe are missing {}", r);
                println!("   do have:");
                for &x in self[sc].iter() {
                    println!("            {}", x);
                }
            }
            assert!(self[sc].contains(&r));
        }
    }

    /// for testing
    pub fn verify_maybe_interacting_includes_everything(&self, r1: Vector3d<Length>) {
        self.verify_subcells_include_everything();
        self.verify_subcells_include_nothing_out_of_place();
        let mi: Vec<_> = self.maybe_interacting_atoms(r1).map(|r| self.put_in_cell(r)).collect();
        let d2s: Vec<_> = self.maybe_interacting_atoms(r1).map(|r| (r-r1).norm2()).collect();
        for &r2 in self.positions.iter() {
            let dist2 = self.closest_distance2(r1, r2);
            if dist2 <= self.well_width*self.well_width {
                if !mi.iter().any(|&rrr| self.closest_distance2(rrr, r2) < 1e-10*units::SIGMA*units::SIGMA) {
                    println!("");
                    println!("  subcell is {} vs {}", self.get_subcell(r2), self.get_subcell(r1));
                    println!("  delta z is {}", r2.z - r1.z);
                    println!("missing [{} {}] {} from:", dist2, d2s.contains(&dist2), r2);
                    for m in mi.iter() {
                        println!("   {}", m);
                    }
                    panic!("We are missing a vector from maybe_interacting_atoms");
                }
            }
        }
    }

    /// for testing
    pub fn verify_maybe_interacting_excluding_includes_everything(&self, r1: Vector3d<Length>, exc: Vector3d<Length>) {
        self.verify_subcells_include_everything();
        self.verify_subcells_include_nothing_out_of_place();
        let mi: Vec<_> = self.maybe_interacting_atoms_excluding(r1, exc).map(|r| self.put_in_cell(r)).collect();
        let d2s: Vec<_> = self.maybe_interacting_atoms_excluding(r1, exc).map(|r| (r-r1).norm2()).collect();
        for &r2 in self.positions.iter().filter(|&&r2| r2 != exc) {
            let dist2 = self.closest_distance2(r1, r2);
            if dist2 <= self.well_width*self.well_width {
                println!("==--> {} from {}", dist2, r2);
                for rr in mi.iter().map(|&r|r).filter(|&rrr| self.closest_distance2(rrr, r2) < 1e-10*units::SIGMA*units::SIGMA) {
                    println!("    ==--> {}", rr);
                }
                if !mi.iter().any(|&rrr| self.closest_distance2(rrr, r2) < 1e-10*units::SIGMA*units::SIGMA) {
                    println!("");
                    println!("==== excluding {} at d2 of {}", exc, self.closest_distance2(exc, r2));
                    println!("==== r1 = {}", r1);
                    println!("==== r2 = {} missing", r2);
                    println!("  box_diagonal is {}", self.box_diagonal);
                    println!("  subcell is {} vs {} / {}", self.get_subcell(r2), self.get_subcell(r1), self.num_subcells);
                    println!("  delta r is {}", r2 - r1);
                    println!("  dist2 = {}", dist2);
                    println!("      mi contains {} with this distance", d2s.iter().filter(|&&d| d == dist2).count());
                    println!("      XX contains {} with this distance",
                             self.positions.iter().filter(|&&r| self.closest_distance2(r1,r) == dist2).count());
                    for x in self[self.get_subcell(r2)].iter() {
                        println!("  subcell has: {}", x);
                    }
                    // for m in mi.iter() {
                    //     if self.closest_distance2(r1, *m) == dist2 {
                    //         println!("**>{}", m);
                    //     } else {
                    //         println!("   {}", m);
                    //     }
                    // }
                    panic!("We are missing a vector from maybe_interacting_atoms_excluding");
                }
            }
        }
    }
}

impl SquareWell {
    fn max_interaction(&self) -> u64 {
        max_balls_within(self.cell.well_width)
    }
    /// Add an atom at a given location.  Returns the change in
    /// energy, or `None` if the atom could not be placed there.
    pub fn add_atom_at(&mut self, r: Vector3d<Length>) -> Option<Energy> {
        let mut e = self.E;
        for r1 in self.cell.maybe_interacting_atoms(r) {
            let dist2 = (r1-r).norm2();
            if dist2 < units::SIGMA*units::SIGMA {
                self.possible_change = Change::None;
                return None;
            } else if dist2 < self.cell.well_width*self.cell.well_width {
                e -= units::EPSILON;
            }
        }
        self.possible_change = Change::Add{ to: r, e };
        Some(e)
    }
    /// Move a specified atom.  Returns the change in energy, or
    /// `None` if the atom could not be placed there.
    pub fn move_atom(&mut self, which: usize, r: Vector3d<Length>) -> Option<Energy> {
        let mut e = self.E;
        let wsqr = self.cell.well_width*self.cell.well_width;
        // let from = self.cell.positions[which];
        let from = unsafe { *self.cell.positions.get_unchecked(which) };
        for r1 in self.cell.maybe_interacting_atoms_excluding(r, from) {
            let dist2 = (r1-r).norm2();
            if dist2 < units::SIGMA*units::SIGMA {
                self.possible_change = Change::None;
                return None;
            }
            if dist2 < wsqr {
                e -= units::EPSILON;
            }
        }
        for r1 in self.cell.maybe_interacting_atoms_excluding(from, from) {
            if (r1-from).norm2() < wsqr {
                e += units::EPSILON;
            }
        }
        self.possible_change = Change::Move{ which, to: r, e };
        Some(e)
    }
    /// Plan to remove the specified atom.  Returns the change in energy.
    pub fn remove_atom_number(&mut self, which: usize) -> Energy {
        let r = self.cell.positions[which];
        let mut e = self.E;
        for r1 in self.cell.maybe_interacting_atoms_excluding(r, r) {
            if self.cell.closest_distance2(r1,r) < self.cell.well_width*self.cell.well_width {
                e += units::EPSILON;
            }
        }
        self.possible_change = Change::Remove{ which, e };
        e
    }
}

impl From<SquareWellParams> for SquareWell {
    fn from(params: SquareWellParams) -> SquareWell {
        SquareWell {
            E: 0.0*units::EPSILON,
            cell: Cell::new(&params._dim, params.well_width*units::SIGMA),
            possible_change: Change::None,
        }
    }
}

impl System for SquareWell {
    fn energy(&self) -> Energy {
        self.E
    }
    fn compute_energy(&self) -> Energy {
        let mut e: Energy = units::EPSILON*0.0;
        for &r1 in self.cell.positions.iter() {
            for r2 in self.cell.maybe_interacting_atoms_excluding(r1, r1) {
                if (r1-r2).norm2() < self.cell.well_width*self.cell.well_width {
                    e -= units::EPSILON;
                }
            }
        }
        e*0.5
    }
    fn delta_energy(&self) -> Option<Energy> {
        Some(units::EPSILON)
    }
    fn greatest_possible_energy(&self) -> Option<Energy> {
        Some(0.0*units::EPSILON)
    }
    fn lowest_possible_energy(&self) -> Option<Energy> {
        Some(-(self.cell.positions.len() as f64)*(self.max_interaction() as f64)*units::EPSILON)
    }
}

impl ConfirmSystem for SquareWell {
    fn confirm(&mut self) {
        match self.possible_change {
            Change::None => (),
            Change::Move{which, to, e} => {
                self.cell.move_atom(which, to);
                self.E = e;
                self.possible_change = Change::None;
            },
            Change::Add{to, e} => {
                self.cell.add_atom_at(to);
                self.E = e;
                self.possible_change = Change::None;
            },
            Change::Remove{which, e} => {
                self.cell.remove_atom(which);
                self.E = e;
                self.possible_change = Change::None;
            },
        }
    }
}

impl GrandSystem for SquareWell {
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
}

impl MovableSystem for SquareWell {
    fn plan_move(&mut self, rng: &mut MyRng, mean_distance: Length) -> Option<Energy> {
        if self.cell.positions.len() > 0 {
            let which = rng.sample(Uniform::new(0, self.cell.positions.len()));
            let to = self.cell.put_in_cell(unsafe { *self.cell.positions.get_unchecked(which) } + rng.vector()*mean_distance);
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
    let mut sw = mk_sw(natoms, 0.3);
    for &r1 in sw.cell.positions.iter() {
        for &r2 in sw.cell.positions.iter() {
            assert_eq!(sw.cell.closest_distance2(r1,r2),
                       sw.cell.unsafe_closest_distance2(r1,r2));
            assert_eq!(sw.cell.closest_distance2(r1,r2),
                       sw.cell.sloppy_closest_distance2(r1,r2));
        }
    }
    let mut rng = MyRng::from_u64(1);
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
    let mut sw = mk_sw(natoms, 0.3);
    let mut rng = MyRng::from_u64(1);
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
    for &r1 in sw.cell.positions.iter() {
        for r2 in sw.cell.maybe_interacting_atoms_excluding(r1, r1) {
            assert_eq!(sw.cell.closest_distance2(r1,r2),
                       sw.cell.unsafe_closest_distance2(r1,r2));
            assert_eq!(sw.cell.closest_distance2(r1,r2),
                       sw.cell.sloppy_closest_distance2(r1,r2));
            assert_eq!(sw.cell.closest_distance2(r1,r2),
                       (r1-r2).norm2());
        }
    }
}

#[test]
fn maybe_interacting_needs_no_shifting_n50() {
    maybe_interacting_needs_no_shifting(50);
}
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
    let mut sw = mk_sw(100, 0.3);
    let mut rng = MyRng::from_u64(1);
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
    let mut sw = mk_sw(natoms, 0.3);
    let mut rng = MyRng::from_u64(1);
    for &r1 in sw.cell.positions.iter() {
        sw.cell.verify_maybe_interacting_excluding_includes_everything(r1, r1);
    }
    for _ in 0..100000 {
        sw.plan_move(&mut rng, Length::new(1.0));
        sw.confirm();
    }
    println!("Finished moving stuff around...");
    for &r1 in sw.cell.positions.iter() {
        sw.cell.verify_maybe_interacting_excluding_includes_everything(r1, r1);
        sw.cell.verify_maybe_interacting_excluding_includes_everything(r1, sw.cell.positions[0]);
        sw.cell.verify_maybe_interacting_excluding_includes_everything(
            r1,
            Vector3d::new(0.5*units::SIGMA, 0.5*units::SIGMA, 0.5*units::SIGMA));
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
fn energy_is_right() {
    let mut sw = SquareWell::from(SquareWellNParams::default());
    assert_eq!(sw.energy(), sw.compute_energy());
    let mut rng = MyRng::from_u64(1);
    for i in 0..1000 {
        println!("making move {}...", i);
        sw.plan_move(&mut rng, Length::new(1.0));
        sw.confirm();
        assert_eq!(sw.energy(), sw.compute_energy());
    }
}

#[cfg(test)]
fn mk_sw(natoms: usize, ff: f64) -> SquareWell {
    let mut param = SquareWellNParams::default();
    param._dim = CellDimensionsGivenNumber::FillingFraction(Unitless::new(ff));
    param.N = natoms;
    SquareWell::from(param)
}
