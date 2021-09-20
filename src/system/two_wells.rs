//! A test model with fake energy.

use super::*;

use dimensioned::Sqrt;
use rand::prelude::*;
use vector3d::Vector3d;
/// The parameters needed to configure a fake model with N dimensions.
///
/// These parameters are normally set via command-line arguments.
#[derive(Serialize, Deserialize, Debug, AutoArgs, Clone, Copy)] //AutoArgs in incompatable with Vec<Length>
#[allow(non_snake_case)]
pub struct Parameters {
    /// the number of dimensions
    pub N: usize,
    /// Ratio of depths of the wells (h_2/h_1)
    pub h2_to_h1: f64,
    /// Energy barrier between wells
    pub barrier_over_h1: f64,
    /// Radius of the second well
    pub r2: Length,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct SystemInvCdf {
    num_points: usize,
    dim: usize,
    r1: Length,
    r2: Length,
    dx0_ball1: Length,
    stencils: Vec<f64>,
}

static NUM_POINTS: usize = 10000;

fn linspace(start: Length, end: Length, numpoints: usize) -> Vec<Length> {
    let last = numpoints as f64 - 1.0;
    let inv_last = 1.0 / last as f64;
    (0..numpoints)
        .map(|i| ((last - i as f64) * start + i as f64 * end) * inv_last)
        .collect()
}

#[allow(non_snake_case)]
impl SystemInvCdf {
    fn new(params: &Parameters) -> Self {
        let dim = params.N;
        let num_points = NUM_POINTS;
        let r1 = Length::new(1.0);
        let r2 = params.r2;

        //  ######## The first axis, ball 1 ########

        let mut stencils = vec![0.; num_points * dim]; // We need a stencil for the first ball
        let numerical_precision_mult = 100;

        // us = np.linspace(-r1, np.sqrt(r1**2 - r2**2) ,(num_points) * numerical_precision_mult)  ## We only need a stencil for the first ball
        let x0: Vec<Length> = linspace(-r1, r1 + 2.0 * r2, num_points);
        let dx0_ball1 = x0[1] - x0[0];

        let dim_power =
            |x_squared: Area| -> f64 { x_squared.sqrt().value_unsafe.powi(dim as i32 - 1) };
        let pdf_x0_nonnormalized = |x0: Length| -> f64 {
            debug_assert!(V(dim - 1).is_finite());
            let pdf = if x0 <= (r1 * r1 - r2 * r2).sqrt() {
                // println!("in first sphere");
                dim_power(r1 * r1 - x0 * x0) * V(dim - 1)
            } else if x0 < r1 + r2 {
                // println!("in cylinder");
                r2.value_unsafe.powi(dim as i32 - 1) * V(dim - 1)
            } else {
                // println!("in hemisphere");
                dim_power(r2 * r2 - (x0 - r1 - r2) * (x0 - r1 - r2)) * V(dim - 1)
            };
            debug_assert!(pdf.is_finite());
            pdf
        };

        let mut val = 0.0;
        for w in 0..num_points - 1 {
            let us = linspace(x0[w], x0[w + 1], numerical_precision_mult);
            let du = us[1] - us[0];
            for i in 0..numerical_precision_mult - 1 {
                // use midpoint method
                debug_assert!(us[i].value_unsafe.is_finite() && us[i + 1].value_unsafe.is_finite());
                val += du.value_unsafe * pdf_x0_nonnormalized(0.5 * (us[i + 1] + us[i]));
            }
            debug_assert!(val.is_finite());
            stencils[w + 1] = val
        }
        assert!(val > 0.);
        for v in stencils.iter_mut() {
            *v /= val;
        }

        // ######## The rest follow a uniform sphere ########

        let xN: Vec<Length> = linspace(Length::new(-1.0), Length::new(1.0), num_points);
        // let dxN = xN[1] - xN[0];
        // println!("dxN {}", dxN);

        let pdf = |x: f64, dim: usize| (1.0 - x * x).powf(0.5 * dim as f64) * V(dim) / V(dim + 1);

        for which_dim in 1..dim {
            let stencil = &mut stencils[which_dim * num_points..(which_dim + 1) * num_points];
            let mut val = 0.0;
            for w in 0..num_points - 1 {
                let us = linspace(xN[w], xN[w + 1], numerical_precision_mult);
                let du = us[1] - us[0];
                for i in 0..numerical_precision_mult - 1 {
                    // use midpoint method
                    debug_assert!(
                        us[i].value_unsafe.is_finite() && us[i + 1].value_unsafe.is_finite()
                    );
                    val += du.value_unsafe
                        * pdf(0.5 * (us[i + 1] + us[i]).value_unsafe, dim - which_dim);
                }
                debug_assert!(val.is_finite());
                stencil[w + 1] = val
            }
            assert!(val > 0.);
            for v in stencil.iter_mut() {
                *v /= val;
            }
        }

        println!("Finished computing inverse cumulative distribution function!");
        SystemInvCdf {
            num_points,
            dim,
            r1,
            r2,
            dx0_ball1,
            stencils,
        }
    }

    fn eval(&self, probability: f64, which_dim: usize) -> Length {
        // L, U, j = find_bin(self.stencils, which_dim*self.num_points, (which_dim+1)*self.num_points-1, probability)
        let stencil =
            &self.stencils[which_dim * self.num_points..(which_dim + 1) * self.num_points];
        let j = binary_search(stencil, probability);
        let L = stencil[j];
        let U = stencil[j + 1];
        if which_dim == 0 {
            let slope = self.dx0_ball1 / (U - L);
            slope * (probability - L) + -self.r1 + j as f64 * self.dx0_ball1
        } else {
            let dx = Length::new(2.0) / (self.num_points as f64 - 1.0);
            let slope = dx / (U - L);
            slope * (probability - L) - Length::new(1.0) + j as f64 * dx
        }
    }

    fn sample(&self, rng: &mut MyRng) -> Vec<Length> {
        let mut sample = Vec::with_capacity(self.dim);
        let x0: Length = self.eval(rng.gen(), 0);
        sample.push(x0);
        let mut R: Length = if x0 <= (self.r1 * self.r1 - self.r2 * self.r2).sqrt() {
            // println!("we're in the first ball");
            (self.r1 * self.r1 - x0 * x0).sqrt()
        } else if x0 < self.r1 + self.r2 {
            // println!("we're in the cylinder");
            self.r2
        } else {
            // println!("we're in the final hemisphere");
            let x2 = x0 - self.r1 - self.r2;
            (self.r2 * self.r2 - x2 * x2).sqrt()
        };
        for i in 2..self.dim {
            let e = self.eval(rng.gen(), i);
            // println!("e: {}", e);
            let x = (R / self.r1) * e;
            sample.push(x);
            R = (R * R - x * x).sqrt();
        }
        sample.push(R * rng.gen_range(-1.0, 1.0));
        sample
    }
}

fn binary_search(values: &[f64], v: f64) -> usize {
    let mut lo = 0;
    let mut hi = values.len() - 1;
    while lo < hi - 1 {
        let mid = (lo + hi) / 2;
        if values[mid] < v {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

// gamma function from https://github.com/rust-lang/rust/issues/18271
#[link(name = "m")]
extern "C" {
    fn tgamma(x: f64) -> f64;
}
fn gamma(x: f64) -> f64 {
    unsafe { tgamma(x) }
}
#[allow(non_snake_case)]
fn V(n: usize) -> f64 {
    return std::f64::consts::PI.powf(0.5 * n as f64) / (gamma(n as f64 * 0.5 + 1.0));
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct Change {
    index: usize,
    /// The new values at index and beyond
    new_x0: Length,
    values: Vector3d<Length>,
}

#[allow(non_snake_case)]
/// A Fake model.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct TwoWells {
    /// The state of the system
    pub position: Vec<Length>,
    /// Squared distance apart from first coordinate
    pub d_squared: Area,
    /// The function itself
    pub parameters: Parameters,
    /// The change we are considering.
    change: Change,
    /// The position of the second well
    well_position: Length,
    /// The inverse cumulative distribution function for computing random positionis
    invcdf: SystemInvCdf,
}

impl From<Parameters> for TwoWells {
    fn from(parameters: Parameters) -> TwoWells {
        // x**2 = b, (w-x)**2/r2**2 = b
        // (w - sqrt(b))**2 = b r2**2
        // w - sqrt(b) = sqrt(b) r2
        // w = sqrt(b)(1 + r2) = sqrt(b)(r1 + r2)
        if parameters.N % 3 != 0 {
            panic!(
                "The number of dimensions {} is not divisible by three!",
                parameters.N
            );
        }
        let well_position = parameters.barrier_over_h1.sqrt() * Length::new(1.)
            + parameters.r2 * (1.0 + parameters.barrier_over_h1 - 1.0 / parameters.h2_to_h1).sqrt();
        let mut position = vec![Length::new(0.0); parameters.N];
        position[0] = Length::new(-0.99);
        let d_squared = position.iter().map(|&x| x * x).sum::<Area>();
        let tw = TwoWells {
            position,
            d_squared,
            change: Change {
                index: 0,
                new_x0: Length::new(0.0),
                values: Vector3d::new(Length::new(0.0), Length::new(0.0), Length::new(0.0)),
            },
            invcdf: SystemInvCdf::new(&parameters),
            well_position,
            parameters,
        };
        for x in linspace(
            Length::new(-1.0),
            Length::new(1.0) + 2.0 * parameters.r2,
            100,
        ) {
            if let Some(energy) = tw.find_energy(x, Area::new(0.0)) {
                println!("{:.4} {:.2}", x, energy);
            } else {
                panic!(
                    "We generated an invalid system?! x0 = {}, d_squared = {}",
                    tw.position[0], tw.d_squared,
                );
            }
        }
        tw
    }
}

impl TwoWells {
    fn find_energy(&self, x0: Length, d_orthog_squared: Area) -> Option<Energy> {
        //r_1 is assumed to be 1 and center_2 is such that the two wells are touching
        let r1 = Length::new(1.0);
        let r2 = self.parameters.r2;
        let rw = self.well_position;

        // let d_orthog_squared = position[1..].iter().map(|&x| x * x).sum::<Area>();

        // let x0 = position[0];
        let x2 = x0 - r1 - r2;
        let xw = x0 - rw;
        let xi = x2 + rw;

        let d_1_squared = d_orthog_squared + x0 * x0;
        let d_2_squared = d_orthog_squared + x2 * x2;
        let d_w_squared = d_orthog_squared + xw * xw;
        let d_i_squared = d_orthog_squared + xi * xi;
        assert!(d_1_squared.value_unsafe.is_finite());

        let e1 = |d2: Area| -> Energy { Energy::new(1.) * (d2 / (r1 * r1) - 1.) };
        let e2 =
            |d2: Area| -> Energy { Energy::new(self.parameters.h2_to_h1) * (d2 / (r2 * r2) - 1.) };

        let e_1 = e1(d_1_squared);
        let e_2 = e2(d_2_squared);
        let e_w = e2(d_w_squared);
        let e_i = e1(d_i_squared);
        // println!("energies are {} {} {} {}", e_1, e_2, e_w, e_i);
        debug_assert!(e_1.value_unsafe.is_finite());
        debug_assert!(e_2.value_unsafe.is_finite());
        debug_assert!(e_w.value_unsafe.is_finite());
        debug_assert!(e_i.value_unsafe.is_finite());

        if d_1_squared <= r1 * r1 {
            // We are in the big sphere, pick whichever is smaller
            Some(if e_1 < e_w { e_1 } else { e_w })
        } else if d_2_squared <= r2 * r2 {
            // We are in the small sphere
            Some(if e_i > e_2 && e_i < Energy::new(0.) {
                e_i
            } else {
                e_2
            })
        } else if d_orthog_squared <= r2 * r2 && x0 > Length::new(0.0) && x0 <= r1 + r2 {
            // We are in the cylinder
            Some(Energy::new(0.0))
        } else {
            None
        }
    }

    fn find_which(&self, x0: Length, d_orthog_squared: Area) -> f64 {
        //r_1 is assumed to be 1 and center_2 is such that the two wells are touching
        let r1 = Length::new(1.0);
        let r2 = self.parameters.r2;
        let rw = self.well_position;

        // let d_orthog_squared = position[1..].iter().map(|&x| x * x).sum::<Area>();

        // let x0 = position[0];
        let x2 = x0 - r1 - r2;
        let xw = x0 - rw;
        let xi = x2 + rw;

        let d_1_squared = d_orthog_squared + x0 * x0;
        let d_2_squared = d_orthog_squared + x2 * x2;
        let d_w_squared = d_orthog_squared + xw * xw;
        let d_i_squared = d_orthog_squared + xi * xi;
        assert!(d_1_squared.value_unsafe.is_finite());

        let e1 = |d2: Area| -> Energy { Energy::new(1.) * (d2 / (r1 * r1) - 1.) };
        let e2 =
            |d2: Area| -> Energy { Energy::new(self.parameters.h2_to_h1) * (d2 / (r2 * r2) - 1.) };

        let e_1 = e1(d_1_squared);
        let e_2 = e2(d_2_squared);
        let e_w = e2(d_w_squared);
        let e_i = e1(d_i_squared);
        // println!("energies are {} {} {} {}", e_1, e_2, e_w, e_i);
        debug_assert!(e_1.value_unsafe.is_finite());
        debug_assert!(e_2.value_unsafe.is_finite());
        debug_assert!(e_w.value_unsafe.is_finite());
        debug_assert!(e_i.value_unsafe.is_finite());

        if d_1_squared <= r1 * r1 {
            // We are in the big sphere, pick whichever is smaller
            if e_1 < e_w {
                0.0 // "big well"
            } else {
                1.0 // "small well"
            }
        } else if d_2_squared <= r2 * r2 {
            // We are in the small sphere
            if e_i > e_2 && e_i < Energy::new(0.) {
                0.0
            } else {
                1.0
            }
        } else if d_orthog_squared <= r2 * r2 && x0 > Length::new(0.0) && x0 < r1 + r2 {
            // We are in the cylinder
            0.0
        } else {
            // We're in an invalid location!
            0.0
        }
    }
}

impl System for TwoWells {
    fn energy(&self) -> Energy {
        self.compute_energy()
    }
    fn compute_energy(&self) -> Energy {
        self.find_energy(
            self.position[0],
            self.d_squared - self.position[0] * self.position[0],
        )
        .expect("position is out of bounds")
    }
    fn randomize(&mut self, rng: &mut MyRng) -> Energy {
        self.position = self.invcdf.sample(rng);
        self.d_squared = self.position.iter().map(|&x| x * x).sum::<Area>();
        if let Some(e) = self.find_energy(
            self.position[0],
            self.d_squared - self.position[0] * self.position[0],
        ) {
            e
        } else {
            println!("We had a roundoff error issue?!");
            self.randomize(rng)
        }
    }
    fn min_moves_to_randomize(&self) -> u64 {
        if self.dimensionality() < 1 {
            1
        } else {
            self.dimensionality() / 3
        }
    }
    fn dimensionality(&self) -> u64 {
        self.position.len() as u64
    }
    fn data_to_collect(&self, _iter: u64) -> Vec<(Interned, f64)> {
        let which_well = Interned::from("which");

        vec![(
            which_well,
            self.find_which(
                self.position[0],
                self.d_squared - self.position[0] * self.position[0],
            ),
        )]
    }
    fn print_debug(&self) {
        print!(
            " x0: {:9.3} {}",
            self.position[0].pretty(),
            if self.find_which(
                self.position[0],
                self.d_squared - self.position[0] * self.position[0]
            ) == 0.0
            {
                "B"
            } else {
                "S"
            }
        );
    }
}

impl ConfirmSystem for TwoWells {
    fn confirm(&mut self) {
        let change = &self.change; // for convenience
        self.d_squared -= self.position[change.index..change.index + 3]
            .iter()
            .map(|&x| x * x)
            .sum::<Area>();
        if change.index != 0 {
            self.d_squared -= self.position[0] * self.position[0];
        }
        self.position[change.index] = change.values.x;
        self.position[change.index + 1] = change.values.y;
        self.position[change.index + 2] = change.values.z;
        self.position[0] = change.new_x0;
        self.d_squared += change.values.norm2();
        if change.index != 0 {
            self.d_squared += change.new_x0 * change.new_x0;
        }
    }
}

impl MovableSystem for TwoWells {
    fn plan_move(&mut self, rng: &mut MyRng, d: Length) -> Option<Energy> {
        use crate::rng::vector;
        use rand::Rng as RandRng;
        let index = 3 * rng.gen_range(0, self.position.len() / 3);
        let old_r = Vector3d::new(
            self.position[index],
            self.position[index + 1],
            self.position[index + 2],
        );
        let r = vector(rng) * d + old_r;
        let d_squared = self.d_squared - old_r.norm2() + r.norm2();
        let x1 = if index == 0 { r.x } else { self.position[0] };
        self.change = Change { index, values: r };
        self.find_energy(x1, d_squared - x1 * x1)
    }
    fn max_size(&self) -> Length {
        Length::new(2.0)
    }
}

#[test]
fn a_test() {
    println!("I am testing");
    assert!(true);
}
