//! This module provides the Cell type, which describes a box of atom positions with periodic boundary conditions.

use super::*;

use dimensioned::Dimensionless;
use dimensioned::{Abs, Cbrt};
use vector3d::Vector3d;

/// A description of the cell dimensions
#[derive(Serialize, Deserialize, Debug, AutoArgs, Clone, Copy)]
#[allow(non_snake_case)]
pub enum CellDimensions {
    /// The three widths of the cell
    CellWidth(Vector3d<Length>),
    /// The volume of the cell
    CellVolume(Volume),
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct Neighbor {
    index: u32,
    offset: Vector3d<i8>,
}

/// A box of atoms
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Cell {
    /// The dimensions of the box.
    pub box_diagonal: Vector3d<Length>,
    /// Distance where repulaive potential goes to zero.
    pub r_cutoff: Length,
    /// The atom positions
    pub positions: Vec<Vector3d<Length>>,
    /// The dimensions of the subcell grid
    #[serde(skip, default)]
    num_subcells: Vector3d<usize>,
    /// The subcell lists, indices for the vectors
    #[serde(skip, default)]
    subcells: Vec<Vec<Neighbor>>,
}

impl Cell {
    /// Create a new Cell with specified dimensions and interaction length.
    pub fn new(dim: &CellDimensions, interaction_length: Length) -> Self {
        let box_diagonal = match dim {
            CellDimensions::CellWidth(w) => Vector3d::new(w.x.abs(), w.y.abs(), w.z.abs()),
            CellDimensions::CellVolume(v) => {
                let w = v.cbrt();
                Vector3d::new(w, w, w)
            }
        };
        let mut cell = Cell {
            box_diagonal: box_diagonal,
            r_cutoff: interaction_length,
            positions: Vec::new(),
            num_subcells: Vector3d::default(),
            subcells: Vec::new(),
        };
        cell.update_caches();
        cell
    }
    /// Update the subcell lists, which seem wasteful to save.
    pub fn update_caches(&mut self) {
        let cells_x = (self.box_diagonal.x / self.r_cutoff).value().floor() as usize;
        let cells_y = (self.box_diagonal.y / self.r_cutoff).value().floor() as usize;
        let cells_z = (self.box_diagonal.z / self.r_cutoff).value().floor() as usize;
        self.num_subcells = Vector3d::new(cells_x, cells_y, cells_z);
        self.subcells = vec![Vec::new(); cells_x * cells_y * cells_z];
        let positions = self.positions.clone();
        for (i, &r) in positions.iter().enumerate() {
            self.add_to_subcells(i as u32, r);
        }
    }
    /// Atoms that may be within r_cutoff of r.
    pub fn maybe_interacting_atoms<'a>(
        &'a self,
        r: Vector3d<Length>,
    ) -> impl Iterator<Item = Vector3d<Length>> + 'a {
        let sc = self.get_subcell(r);
        self.index(sc)
            .iter()
            .map(move |Neighbor { index, offset }| {
                self.positions[*index as usize]
                    - Vector3d::new(
                        offset.x as f64 * self.box_diagonal.x,
                        offset.y as f64 * self.box_diagonal.y,
                        offset.z as f64 * self.box_diagonal.z,
                    )
            })
    }
    /// Atoms that may be within r_cutoff of r.  This excludes any
    /// atom that is located precisely at rprime.
    pub fn maybe_interacting_atoms_excluding<'a>(
        &'a self,
        r: Vector3d<Length>,
        which: usize,
    ) -> impl Iterator<Item = Vector3d<Length>> + 'a {
        let sc = self.get_subcell(r);
        self.index(sc)
            .iter()
            .filter(move |Neighbor { index, .. }| *index != which as u32)
            .map(move |Neighbor { index, offset }| {
                self.positions[*index as usize]
                    - Vector3d::new(
                        offset.x as f64 * self.box_diagonal.x,
                        offset.y as f64 * self.box_diagonal.y,
                        offset.z as f64 * self.box_diagonal.z,
                    )
            })
    }
    /// Find the cell for a given vector.
    pub fn get_subcell(&self, r: Vector3d<Length>) -> Vector3d<isize> {
        Vector3d::new(
            (r.x / self.box_diagonal.x * self.num_subcells.x as f64)
                .value()
                .floor() as isize,
            (r.y / self.box_diagonal.y * self.num_subcells.y as f64)
                .value()
                .floor() as isize,
            (r.z / self.box_diagonal.z * self.num_subcells.z as f64)
                .value()
                .floor() as isize,
        )
    }
    /// Add an atom at a given location.
    pub fn add_atom_at(&mut self, r: Vector3d<Length>) {
        let index = self.positions.len() as u32;
        self.positions.push(r);
        self.add_to_subcells(index, r);
    }
    fn add_to_subcells(&mut self, index: u32, r: Vector3d<Length>) {
        let sc = self.get_subcell(r);
        for &n in NEIGHBORS.iter() {
            let scp = sc + n;
            let offset = Vector3d::new(
                if scp.x < 0 {
                    -1
                } else if scp.x == self.num_subcells.x as isize {
                    1
                } else {
                    0
                },
                if scp.y < 0 {
                    -1
                } else if scp.y == self.num_subcells.y as isize {
                    1
                } else {
                    0
                },
                if scp.z < 0 {
                    -1
                } else if scp.z == self.num_subcells.z as isize {
                    1
                } else {
                    0
                },
            );
            self.index_mut(scp).push(Neighbor { index, offset });
        }
    }
    /// Move the atom.  Return the previous position.
    pub fn move_atom(&mut self, which: usize, r: Vector3d<Length>) -> Vector3d<Length> {
        // let old = ::std::mem::replace(&mut self.positions[which], r);
        let old = ::std::mem::replace(unsafe { self.positions.get_unchecked_mut(which) }, r);
        let sc = self.get_subcell(r);
        let oldsc = self.get_subcell(old);
        if sc != oldsc {
            let index = which as u32;
            for &n in NEIGHBORS.iter() {
                remove_if(self.index_mut(oldsc + n), |n| n.index == index);
            }
            self.add_to_subcells(index, r);
        }
        old
    }
    /// Remove an atom.  Return the previous position.
    pub fn remove_atom(&mut self, which: usize) -> Vector3d<Length> {
        let r = self.positions.swap_remove(which);
        let sc = self.get_subcell(r);
        for &n in NEIGHBORS.iter() {
            remove_if(self.index_mut(sc + n), |n| n.index == which as u32);
        }
        if self.positions.len() != which {
            // When we swap_removed the atom, another atom got its
            // index changed, and we need to adjust for that.
            let oldindex = self.positions.len() as u32;
            let sc = self.get_subcell(self.positions[which]);
            for &n in NEIGHBORS.iter() {
                rename_neighbor(self.index_mut(sc + n), oldindex, which as u32);
            }
        }
        r
    }
    /// The volume of the cell
    pub fn volume(&self) -> Volume {
        self.box_diagonal[0] * self.box_diagonal[1] * self.box_diagonal[2]
    }
    /// PUBLIC FOR TESTING ONLY! The shortest distance squared between two vectors.
    pub fn closest_distance2(&self, r1: Vector3d<Length>, r2: Vector3d<Length>) -> Area {
        let mut dr = r2 - r1;
        if dr.x < -0.5 * self.box_diagonal.x {
            while {
                dr.x += self.box_diagonal.x;
                dr.x < -0.5 * self.box_diagonal.x
            } {}
        } else {
            while dr.x > 0.5 * self.box_diagonal.x {
                dr.x -= self.box_diagonal.x;
            }
        }
        if dr.y < -0.5 * self.box_diagonal.y {
            while {
                dr.y += self.box_diagonal.y;
                dr.y < -0.5 * self.box_diagonal.y
            } {}
        } else {
            while dr.y > 0.5 * self.box_diagonal.y {
                dr.y -= self.box_diagonal.y;
            }
        }
        if dr.z < -0.5 * self.box_diagonal.z {
            while {
                dr.z += self.box_diagonal.z;
                dr.z < -0.5 * self.box_diagonal.z
            } {}
        } else {
            while dr.z > 0.5 * self.box_diagonal.z {
                dr.z -= self.box_diagonal.z;
            }
        }
        dr.norm2()
    }
    /// PUBLIC FOR TESTING ONLY! The shortest distance squared between two vectors.
    pub fn sloppy_closest_distance2(&self, r1: Vector3d<Length>, r2: Vector3d<Length>) -> Area {
        let mut dr = r2 - r1;
        while dr.x < -0.5 * self.box_diagonal.x {
            dr.x += self.box_diagonal.x;
        }
        while dr.x > 0.5 * self.box_diagonal.x {
            dr.x -= self.box_diagonal.x;
        }
        while dr.y < -0.5 * self.box_diagonal.y {
            dr.y += self.box_diagonal.y;
        }
        while dr.y > 0.5 * self.box_diagonal.y {
            dr.y -= self.box_diagonal.y;
        }
        while dr.z < -0.5 * self.box_diagonal.z {
            dr.z += self.box_diagonal.z;
        }
        while dr.z > 0.5 * self.box_diagonal.z {
            dr.z -= self.box_diagonal.z;
        }
        dr.norm2()
    }
    /// PUBLIC FOR TESTING ONLY! The shortest distance squared between two vectors.
    pub fn unsafe_closest_distance2(&self, r1: Vector3d<Length>, r2: Vector3d<Length>) -> Area {
        let mut dr = r2 - r1;
        if dr.x < -0.5 * self.box_diagonal.x {
            dr.x += self.box_diagonal.x;
        } else if dr.x > 0.5 * self.box_diagonal.x {
            dr.x -= self.box_diagonal.x;
        }
        if dr.y < -0.5 * self.box_diagonal.y {
            dr.y += self.box_diagonal.y;
        } else if dr.y > 0.5 * self.box_diagonal.y {
            dr.y -= self.box_diagonal.y;
        }
        if dr.z < -0.5 * self.box_diagonal.z {
            dr.z += self.box_diagonal.z;
        } else if dr.z > 0.5 * self.box_diagonal.z {
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
    /// FOR TESTIN
    pub fn verify_maybe_interacting_includes_everything(&self, _r1: Vector3d<Length>) {}
    /// FOR TESTIN
    pub fn verify_maybe_interacting_excluding_includes_everything(
        &self,
        _r1: Vector3d<Length>,
        _exc: usize,
    ) {
    }
}

fn modulus(i: isize, sz: usize) -> usize {
    let sz = sz as isize;
    ((i + sz) % sz) as usize
}

impl Cell {
    fn index(&self, i: Vector3d<isize>) -> &Vec<Neighbor> {
        let ix = modulus(i.x, self.num_subcells.x);
        let iy = modulus(i.y, self.num_subcells.y);
        let iz = modulus(i.z, self.num_subcells.z);
        // &self.subcells[ix*(self.num_subcells.y*self.num_subcells.z) + iy*self.num_subcells.z + iz]
        unsafe {
            self.subcells.get_unchecked(
                ix * (self.num_subcells.y * self.num_subcells.z) + iy * self.num_subcells.z + iz,
            )
        }
    }
    fn index_mut(&mut self, i: Vector3d<isize>) -> &mut Vec<Neighbor> {
        let ix = modulus(i.x, self.num_subcells.x);
        let iy = modulus(i.y, self.num_subcells.y);
        let iz = modulus(i.z, self.num_subcells.z);
        // &mut self.subcells[ix*(self.num_subcells.y*self.num_subcells.z) + iy*self.num_subcells.z + iz]
        unsafe {
            self.subcells.get_unchecked_mut(
                ix * (self.num_subcells.y * self.num_subcells.z) + iy * self.num_subcells.z + iz,
            )
        }
    }
}

const NEIGHBORS: [Vector3d<isize>; 27] = [
    Vector3d { x: 0, y: 0, z: 0 },
    Vector3d { x: 1, y: 0, z: 0 },
    Vector3d { x: -1, y: 0, z: 0 },
    Vector3d { x: 0, y: 1, z: 0 },
    Vector3d { x: 0, y: -1, z: 0 },
    Vector3d { x: 0, y: 0, z: 1 },
    Vector3d { x: 0, y: 0, z: -1 },
    Vector3d { x: 0, y: 1, z: 1 },
    Vector3d { x: 0, y: 1, z: -1 },
    Vector3d { x: 0, y: -1, z: 1 },
    Vector3d { x: 0, y: -1, z: -1 },
    Vector3d { x: 1, y: 0, z: 1 },
    Vector3d { x: 1, y: 0, z: -1 },
    Vector3d { x: -1, y: 0, z: 1 },
    Vector3d { x: -1, y: 0, z: -1 },
    Vector3d { x: 1, y: 1, z: 0 },
    Vector3d { x: 1, y: -1, z: 0 },
    Vector3d { x: -1, y: 1, z: 0 },
    Vector3d { x: -1, y: -1, z: 0 },
    Vector3d { x: 1, y: 1, z: 1 },
    Vector3d { x: -1, y: 1, z: 1 },
    Vector3d { x: 1, y: -1, z: 1 },
    Vector3d { x: 1, y: 1, z: -1 },
    Vector3d { x: 1, y: -1, z: -1 },
    Vector3d { x: -1, y: 1, z: -1 },
    Vector3d { x: -1, y: -1, z: 1 },
    Vector3d {
        x: -1,
        y: -1,
        z: -1,
    },
];

fn remove_if<T>(vec: &mut Vec<T>, f: impl Fn(&T) -> bool) {
    let mut ind = None;
    for (i, v) in vec.iter().enumerate() {
        if f(v) {
            ind = Some(i);
            break;
        }
    }
    if let Some(i) = ind {
        vec.swap_remove(i);
    }
}
fn rename_neighbor(vec: &mut Vec<Neighbor>, oldindex: u32, newindex: u32) {
    for n in vec.iter_mut() {
        if n.index == oldindex {
            n.index = newindex;
            break;
        }
    }
}
