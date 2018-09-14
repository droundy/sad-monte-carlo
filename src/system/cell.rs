//! This module provides the Cell type, which describes a box of atom positions.

use super::*;

use dimensioned::Dimensionless;
use dimensioned::{Cbrt,Abs};
use vector3d::Vector3d;

/// A description of the cell dimensions
#[derive(Serialize, Deserialize, Debug, ClapMe)]
#[allow(non_snake_case)]
pub enum CellDimensions {
    /// The three widths of the cell
    CellWidth(Vector3d<Length>),
    /// The volume of the cell
    CellVolume(Volume),
}

/// A box of atoms
#[derive(Serialize, Deserialize, Debug)]
pub struct Cell {
    /// The dimensions of the box.
    pub box_diagonal: Vector3d<Length>,
    /// The well width.
    pub well_width: Length,
    /// The atom positions
    pub positions: Vec<Vector3d<Length>>,
    /// The dimensions of the subcell grid
    num_subcells: Vector3d<usize>,
    /// The subcell lists
    subcells: Vec<Vec<Vector3d<Length>>>
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
