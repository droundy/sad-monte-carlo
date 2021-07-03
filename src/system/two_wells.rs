//! A test model with fake energy.

use super::*;

use dimensioned::Sqrt;
use rand::prelude::*;
/// The parameters needed to configure a fake model with N dimensions.
///
/// These parameters are normally set via command-line arguments.
#[derive(Serialize, Deserialize, Debug, AutoArgs, Clone)] //AutoArgs in incompatable with Vec<Length>
#[allow(non_snake_case)]
pub struct Parameters {
    /// the number of dimensions
    pub N: usize,
    /// Ratio of depths of the wells (h_2/h_1)
    pub h2_to_h1: f64,
    /// Radius of the second well
    pub r2: Length,
}

struct SystemInvCdf {
    num_points: usize,
    dim: usize,
    r1: Length,
    r2: Length,
    dx1_ball1: Length,
    stencils: Vec<f64>,
}

static NUM_POINTS: usize = 10000;

impl SystemInvCdf {
    fn new(params: &Parameters) -> Self {
        let dim = params.N;
        let num_points = NUM_POINTS;
        let r1 = Length::new(1.0);
        let r2 = params.r2;
        let dx1_ball1 = (r1 + (r1 * r1 - r2 * r2).sqrt()) / num_points as f64;

        //  ######## The first axis, ball 1 ########

        let mut stencils = vec![0.; num_points * dim]; // We need a stencil for the first ball
        let numerical_precision_mult = 100;

        // us = np.linspace(-r1, np.sqrt(r1**2 - r2**2) ,(num_points) * numerical_precision_mult)  ## We only need a stencil for the first ball
        let mut us: Vec<Length> = vec![Length::new(0.0); num_points * numerical_precision_mult];
        for i in 0..us.len() {
            us[i] =
                -r1 + i as f64 * ((r1 * r1 - r2 * r2).sqrt() + r1) * (1.0 / (us.len() as f64 - 1.0))
        }
        let du = us[1] - us[0];

        let dim_power =
            |x_squared: Area| -> f64 { x_squared.sqrt().value_unsafe.powi(dim as i32 - 1) };
        let pdf_x1_nonnormalized = |x1: Length| -> f64 {
            if x1 <= (r1 * r1 - r2 * r2).sqrt() {
                dim_power(r1 * r1 - x1 * x1) * V(dim - 1)
            } else if x1 < r1 + r2 {
                r2.value_unsafe.powi(dim as i32 - 1) * V(dim - 1)
            } else {
                dim_power(r2 * r2 - (x1 - r1 - r2) * (x1 - r1 - r2)) * V(dim - 1)
            }
        };

        let mut val = 0.0;
        let mut w = 1;
        for i in 1..(num_points - 1) * numerical_precision_mult {
            if i % numerical_precision_mult == 0 {
                stencils[w] = val;
                w += 1;
            }
            val += (pdf_x1_nonnormalized(us[i - 1])
                + 4. * pdf_x1_nonnormalized((us[i] + us[i - 1]) * 0.5)
                + pdf_x1_nonnormalized(us[i]))
                * du.value_unsafe
                / 6.0;
        }
        stencils[w] = val;
        let V_ball1 = val;
        for v in stencils.iter_mut() {
            *v /= val;
        }

        // ######## Find the volume of the cylinder, second ball, and total volume ########

        let V_ball2 = V(dim) * r2.value_unsafe.powi(dim as i32) * 0.5;
        let V_cyl = (r1 + 2. * r2 - (r1 * r1 - r2 * r2).sqrt()).value_unsafe
            * V(dim - 1)
            * r2.value_unsafe.powi(dim as i32);
        let V_tot = V_ball1 + V_ball2 + V_cyl;

        // ######## calculate Probabilities of being in Ball1, the cylinder, or Ball 2 ########

        let P_ball1 = V_ball1 / V_tot;
        let P_cyl = V_cyl / V_tot;
        let P_ball2 = V_ball2 / V_tot;

        // ######## The rest follow a uniform sphere ########

        //         us = np.linspace(-1,1,(num_points-1) * numerical_precision_mult)
        let mut us: Vec<f64> = vec![0.0; (num_points - 1) * numerical_precision_mult];
        for i in 0..us.len() {
            us[i] = -1.0 + i as f64 * (1.0 / (us.len() as f64 - 1.0));
        }
        let du = us[1] - us[0];

        let pdf = |x: f64, dim: usize| (1.0 - x * x).powf(0.5 * dim as f64) * V(dim) / V(dim + 1);

        for N in 1..dim {
            val = 0.0;
            w = 1;
            for i in 1..(num_points - 1) * numerical_precision_mult {
                if i % numerical_precision_mult == 0 {
                    stencils[N * num_points + w] = val;
                    w += 1;
                }
                val += (pdf(us[i - 1], dim - N)
                    + 4. * pdf((us[i] + us[i - 1]) * 0.5, dim - N)
                    + pdf(us[i], dim - N))
                    * du
                    / 6.;
            }
            stencils[N * num_points + w] = 1.; // val
        }

        SystemInvCdf {
            num_points,
            dim,
            r1,
            r2,
            dx1_ball1,
            stencils,
        }
    }

    fn eval(&self, probability: f64, which_dim: usize) -> Length {
        // L, U, j = find_bin(self.stencils, which_dim*self.num_points, (which_dim+1)*self.num_points-1, probability)
        let stencil =
            &self.stencils[which_dim * self.num_points..(which_dim + 1) * self.num_points];
        let j = binary_search(stencil, probability);
        let L = stencil[j];
        let U = stencil[j+1];
        if which_dim == 0 {
            let slope = self.dx1_ball1/(U - L);
            slope * (probability - L) + -self.r1 + j as f64*self.dx1_ball1
        } else {
            let dx = Length::new(2.0/self.num_points as f64);
            let slope = dx/(U-L);
            slope * (probability - L) - Length::new(1.0) + j as f64*dx
        }
    }
}

fn binary_search(values: &[f64], v: f64) -> usize {
    let mut lo = 0;
    let mut hi = values.len() - 1;
    while lo < hi - 1 {
        let mid = (lo + hi)/2;
        if values[mid] < v {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

//     def print_data(self, which_dim):
//         print(self.stencils[which_dim * self.num_points : (which_dim+1) * self.num_points])

//     def sample(self):
//         sample = np.zeros(self.dim)
//         ### Determine if in ball 1 ###
//         random_num = np.random.uniform(0,1)

//         if(random_num < self.P_ball1):
//             sample[0] = self.eval(np.random.uniform(0,1), 0)
//             R = np.sqrt(self.r1**2 - sample[0]**2)
//         ### Determine if in the cylinder ###
//         elif(random_num < self.P_ball1 + self.P_cyl):
//             sample[0] = np.random.uniform(np.sqrt(self.r1**2 - self.r2**2), self.r1 + self.r2)
//             R = self.r2
//         ### Determine if in ball 2 ###
//         else:
//             sample[0] = np.abs( self.r2 * self.eval(np.random.uniform(0,1), 1) )
//             R = np.sqrt(self.r2**2 - sample[0]**2)
//             sample[0] += self.r1 + self.r2

//         for i in range(2,self.dim):
//             sample[i-1] = self.eval(np.random.uniform(0,1), i)
//             sample[i-1] *= R
//             R = np.sqrt(R**2 - sample[i-1]**2)

//         sample[self.dim - 1] = R * np.random.uniform(-1,1)

//         return sample

// gamma function from https://github.com/rust-lang/rust/issues/18271
#[link(name = "m")]
extern "C" {
    fn tgamma(x: f64) -> f64;
}
fn gamma(x: f64) -> f64 {
    unsafe { tgamma(x) }
}
fn V(n: usize) -> f64 {
    return std::f64::consts::PI.powf(0.5 * n as f64) / (gamma(n as f64 * 0.5 + 1.0));
}

#[allow(non_snake_case)]
/// A Fake model.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct TwoWells {
    /// The state of the system
    pub position: Vec<Length>,
    /// The function itself
    pub parameters: Parameters,
    /// The last change we made (and might want to undo).
    possible_change: Vec<Length>,
}

impl From<Parameters> for TwoWells {
    fn from(parameters: Parameters) -> TwoWells {
        TwoWells {
            position: vec![Length::new(0.5); parameters.N],
            possible_change: Vec::new(),
            parameters,
        }
    }
}

impl TwoWells {
    fn find_energy(&self, position: &[Length]) -> Energy {
        //r_1 is assumed to be 1 and center_2 is such that the two wells are touching
        let r1 = Length::new(1.0);

        let mut d_1_squared = Area::new(0.);
        let mut d_2_squared = Area::new(0.);

        let center_2 = (Length::new(1.) + self.parameters.r2) / ((self.parameters.N as f64).sqrt());

        for i in 0..self.parameters.N {
            d_1_squared += position[i] * position[i];
            d_2_squared += (position[i] - center_2) * (position[i] - center_2);
        }

        let e_1 = Energy::new(1.) * (d_1_squared - r1 * r1) / (r1 * r1);
        let e_2 = Energy::new(self.parameters.h2_to_h1)
            * (d_2_squared / (self.parameters.r2 * self.parameters.r2) - 1.);

        if e_1 < e_2 {
            e_1
        } else {
            e_2
        }
    }
}

impl System for TwoWells {
    fn energy(&self) -> Energy {
        self.find_energy(&self.position)
    }
    fn compute_energy(&self) -> Energy {
        self.energy()
    }
    fn randomize(&mut self, rng: &mut MyRng) -> Energy {
        for x in self.position.iter_mut() {
            *x = Length::new(rng.gen_range(-1.0, 1.0));
        }
        self.energy()
    }
    fn min_moves_to_randomize(&self) -> u64 {
        self.dimensionality() // FIXME /3
    }
    fn dimensionality(&self) -> u64 {
        self.position.len() as u64
    }
}

impl ConfirmSystem for TwoWells {
    fn confirm(&mut self) {
        self.position = self.possible_change.clone();
    }
}

impl MovableSystem for TwoWells {
    fn plan_move(&mut self, rng: &mut MyRng, d: Length) -> Option<Energy> {
        // FIXME move three coordinates at once.
        let i = rng.gen_range(0, self.position.len());
        self.possible_change = self.position.clone();
        let v: f64 = rng.sample(rand_distr::StandardNormal);
        self.possible_change[i] += v * d;
        if self.possible_change[i] >= Length::new(1.0)
            || self.possible_change[i] <= Length::new(-1.0)
        {
            return None;
        }
        Some(self.find_energy(&self.possible_change))
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
