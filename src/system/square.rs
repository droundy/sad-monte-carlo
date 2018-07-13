//! A square well fluid.

use super::*;

use dimensioned::Dimensionless;
use vector3d::Vector3d;
use clapme::ClapMe;
use clapme;

#[derive(Serialize, Deserialize, Debug, ClapMe)]
pub struct SquareWellParams {
    well_width: Unitless,
    box_diagonal: Vector3d<Length>,
}

#[allow(non_snake_case)]
/// A square well fluid.
#[derive(Serialize, Deserialize, Debug)]
pub struct SquareWell {
    /// The dimensionless well width.
    well_width: Length,
    /// The atom positions
    positions: Vec<Vector3d<Length>>,
    /// The energy of the system
    E: Energy,
    /// The dimensions of the box.
    box_diagonal: Vector3d<Length>,
    /// The last change we made (and might want to undo).
    last_change: Change,
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
enum Change {
    Move { which: usize, from: Vector3d<Length>, de: Energy },
    Add { which: usize, de: Energy },
    Remove { from: Vector3d<Length>, de: Energy },
    None,
}

impl SquareWell {
    /// Create a new square well fluid.
    pub fn new(params: SquareWellParams) -> SquareWell {
        let mut box_diagonal = params.box_diagonal;
        if box_diagonal.x < units::SIGMA*0.0 {
            box_diagonal.x = -box_diagonal.x;
        }
        if box_diagonal.y < units::SIGMA*0.0 {
            box_diagonal.y = -box_diagonal.y;
        }
        if box_diagonal.z < units::SIGMA*0.0 {
            box_diagonal.z = -box_diagonal.z;
        }
        SquareWell {
            well_width: params.well_width*units::SIGMA,
            positions: Vec::new(),
            E: 0.0*units::EPSILON,
            box_diagonal: box_diagonal,
            last_change: Change::None,
        }
    }
    fn max_interaction(&self) -> u64 {
        max_balls_within(self.well_width)
    }
    /// Add an atom at a given location.  Returns the change in
    /// energy, or `None` if the atom could not be placed there.
    pub fn add_atom_at(&mut self, r: Vector3d<Length>) -> Option<Energy> {
        let mut de = units::EPSILON*0.0;
        for &r1 in self.positions.iter() {
            let dist2 = self.closest_distance2(r1,r);
            if dist2 < units::SIGMA*units::SIGMA {
                self.last_change = Change::None;
                return None;
            } else if dist2 < self.well_width*self.well_width {
                de -= units::EPSILON;
            }
        }
        self.positions.push(r);
        self.last_change = Change::Add{ which: self.positions.len()-1, de: de };
        Some(de)
    }
    /// Remove the specified atom.  Returns the change in energy.
    pub fn remove_atom_number(&mut self, which: usize) -> Energy {
        let r = self.positions.swap_remove(which);
        let mut de = units::EPSILON*0.0;
        for &r1 in self.positions.iter() {
            if self.closest_distance2(r1,r) < self.well_width*self.well_width {
                de += units::EPSILON;
            }
        }
        self.last_change = Change::Remove{ from: r, de: de };
        de
    }
    fn closest_distance2(&self, r1: Vector3d<Length>, r2: Vector3d<Length>) -> Area {
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
}

impl System for SquareWell {
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

impl UndoSystem for SquareWell {
    fn undo(&mut self) {
        match self.last_change {
            Change::None => (),
            Change::Move{which, from, de} => {
                let old = self.positions[which];
                self.positions[which] = from;
                self.E -= de;
                self.last_change = Change::Move {
                    which: which,
                    from: old,
                    de: -de,
                };
            },
            Change::Add{which, de} => {
                let old = self.positions.swap_remove(which);
                self.E -= de;
                self.last_change = Change::Remove {
                    from: old,
                    de: -de,
                };
            },
            Change::Remove{from, de} => {
                self.positions.push(from);
                self.E -= de;
                self.last_change = Change::Add {
                    which: self.positions.len()-1,
                    de: -de,
                };
            },
        }
    }
}

impl GrandSystem for SquareWell {
    fn add_atom(&mut self) -> Option<Energy> {
        None
    }
    fn remove_atom(&mut self) -> Energy {
        Energy::new(0.0)
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
