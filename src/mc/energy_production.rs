//! An implementation of SAD (Statistical Association, Dynamical version).
//!
//! This version using a [Binning]() approach so that we can use
//! something nicer than a histogram.  The idea is also to use this
//! for a 2D histogram, which will be even fancier.

#![allow(non_snake_case)]

use super::*;
use crate::system::*;

use rand::{Rng, SeedableRng};
use std::default::Default;

/// The parameters needed to configure a simulation.
#[derive(Debug, AutoArgs, Clone)]
pub struct MCParams {
    /// The base filename for which we are running the production
    pub base: Option<String>,
    /// Subdivide each bin by this much.
    pub subdivide: Option<usize>,
    /// The seed for the random number generator.
    pub seed: Option<u64>,
    /// report input
    pub _movies: plugin::MovieParams,
    /// report input
    pub _save: plugin::SaveParams,
    /// report input
    pub _report: plugin::ReportParams,
}

impl Default for MCParams {
    fn default() -> Self {
        MCParams {
            base: None,
            seed: None,
            subdivide: None,
            _save: plugin::SaveParams::default(),
            _movies: plugin::MovieParams::default(),
            _report: plugin::ReportParams::default(),
        }
    }
}

/// A simulation with many replicas
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct MC<S> {
    /// The random number generator.
    pub rng: crate::rng::MyRng,
    /// Where to save the resume file.
    pub save_as: std::path::PathBuf,
    /// The number of MC moves we have done
    pub moves: u64,
    /// Number of rejected moves
    pub rejected_moves: u64,

    /// The system itself
    pub system: S,

    /// The energy boundaries
    pub energy_boundaries: Vec<Energy>,
    /// The lnw
    pub lnw: Vec<f64>,

    /// The histogram
    pub histogram: Vec<u64>,
    /// The total energy in each bin
    pub total_energy: Vec<Energy>,

    /// How frequently to save...
    save: plugin::Save,
    /// Movie state
    movie: plugin::Movie,
    /// Movie state
    report: plugin::Report,
}

fn read_f64(fname: String) -> Result<Vec<f64>, std::io::Error> {
    use std::io::Read;
    let mut f = std::fs::File::open(&fname)?;

    let mut v = Vec::new();
    f.read_to_end(&mut v)?;

    Ok(v.split(|c| b"\n ".contains(c))
        .map(|s| String::from_utf8_lossy(s))
        .flat_map(|s| s.parse().ok())
        .collect())
}

fn divide_boundaries(bounds: &[f64], subdivide: usize) -> Vec<f64> {
    let mut eb = Vec::new();
    let de = (bounds[0] - bounds[1]) / subdivide as f64;
    for i in (1..subdivide).rev() {
        eb.push(bounds[0] + i as f64 * de);
    }
    for b in 0..bounds.len() - 1 {
        let de = (bounds[b] - bounds[b + 1]) / subdivide as f64;
        for i in 0..subdivide {
            eb.push(bounds[b] - i as f64 * de);
        }
    }
    let de = (bounds[bounds.len() - 2] - bounds[bounds.len() - 1]) / subdivide as f64;
    for i in 0..subdivide {
        eb.push(bounds[bounds.len() - 1] - i as f64 * de);
    }
    eb
}

#[test]
fn bounds_subdivision() {
    let bounds = &[1.0, 0.0];
    let half = divide_boundaries(bounds, 2);
    let quarter = divide_boundaries(bounds, 4);
    assert_eq!(&half, &[1.5, 1.0, 0.5, 0.0, -0.5]);
    assert_eq!(
        &quarter,
        &[1.75, 1.5, 1.25, 1.0, 0.75, 0.5, 0.25, 0.0, -0.25, -0.5, -0.75]
    );
}

impl<S: Clone + ConfirmSystem + MovableSystem + serde::Serialize + serde::de::DeserializeOwned>
    MC<S>
{
    fn from_params(params: MCParams, system: S, save_as: std::path::PathBuf) -> Self {
        let rng = crate::rng::MyRng::seed_from_u64(params.seed.unwrap_or(0));

        let base = if let Some(b) = params.base {
            b
        } else {
            save_as.with_extension("").to_str().unwrap().to_string()
        };
        let mut energy_boundaries = read_f64(format!("{}-energy-boundaries.dat", base)).unwrap();
        let mut lnw = read_f64(format!("{}-lnw.dat", base)).unwrap();
        if let Some(subdivide) = params.subdivide {
            lnw = lnw
                .iter()
                .cloned()
                .flat_map(|v| vec![v; subdivide])
                .collect();
            energy_boundaries = divide_boundaries(&energy_boundaries, subdivide);
            assert_eq!(lnw.len(), energy_boundaries.len() + 1);
        }
        let energy_boundaries: Vec<Energy> = energy_boundaries
            .iter()
            .cloned()
            .map(|x| Energy::new(x))
            .collect();
        MC {
            moves: 0,
            rejected_moves: 0,

            system,
            rng,
            histogram: vec![0; lnw.len()],
            total_energy: vec![Energy::new(0.0); lnw.len()],
            energy_boundaries,
            lnw,
            save_as: save_as,
            save: plugin::Save::from(params._save),
            movie: plugin::Movie::from(params._movies),
            report: plugin::Report::from(params._report),
        }
    }

    /// Create a new simulation from command-line flags.
    pub fn from_args<A: AutoArgs + Into<S>>() -> Self {
        println!("git version: {}", VERSION);
        match <Params<MCParams, A>>::from_args() {
            Params::_Params {
                _sys,
                _mc,
                save_as,
                num_threads,
            } => {
                if let Some(num_threads) = num_threads {
                    rayon::ThreadPoolBuilder::new()
                        .num_threads(num_threads)
                        .build_global()
                        .unwrap()
                }
                if let Some(ref save_as) = save_as {
                    if let Ok(f) = ::std::fs::File::open(save_as) {
                        let mut s = match save_as.extension().and_then(|x| x.to_str()) {
                            Some("yaml") => serde_yaml::from_reader::<_, Self>(&f)
                                .expect("error parsing save-as file"),
                            Some("json") => serde_json::from_reader::<_, Self>(&f)
                                .expect("error parsing save-as file"),
                            Some("cbor") => serde_cbor::from_reader::<Self, _>(&f)
                                .expect("error parsing save-as file"),
                            _ => panic!("I don't know how to read file {:?}", f),
                        };
                        println!("Resuming from file {:?}", save_as);
                        s.system.update_caches();
                        s.report.update_from(_mc._report);
                        s.save.update_from(_mc._save);
                        return s;
                    } else {
                        return Self::from_params(_mc, _sys.into(), save_as.clone());
                    }
                }
                let save_as = save_as.unwrap_or(::std::path::PathBuf::from("resume.yaml"));
                Self::from_params(_mc, _sys.into(), save_as)
            }
            Params::ResumeFrom(p) => {
                let f = ::std::fs::File::open(&p).expect(&format!("error reading file {:?}", &p));
                match p.extension().and_then(|x| x.to_str()) {
                    Some("yaml") => {
                        serde_yaml::from_reader(&f).expect("error reading checkpoint?!")
                    }
                    Some("json") => {
                        serde_json::from_reader(&f).expect("error reading checkpoint?!")
                    }
                    Some("cbor") => {
                        serde_cbor::from_reader(&f).expect("error reading checkpoint?!")
                    }
                    _ => panic!("I don't know how to read file {:?}", f),
                }
            }
        }
    }
    /// Create a simulation checkpoint.
    pub fn checkpoint(&mut self) {
        self.report.print(self.moves, self.moves);
        println!(
            "    Current energy {:.5} rejected {:.2}%",
            self.system.energy().pretty(),
            crate::prettyfloat::PrettyFloat(self.rejected_moves as f64 / self.moves as f64 * 100.0)
        );
        let unexplored = self.histogram.iter().cloned().filter(|&h| h == 0).count();
        if unexplored > 0 {
            println!(
                "    {:.4} bins remain unexplored",
                crate::prettyfloat::PrettyFloat(unexplored as f64)
            );
        } else {
            let minimum = self.histogram.iter().cloned().min().unwrap();
            println!(
                "    saw each bin at least {:.2} times",
                crate::prettyfloat::PrettyFloat(minimum as f64)
            );
        }

        let f = AtomicFile::create(&self.save_as)
            .expect(&format!("error creating file {:?}", self.save_as));
        match self.save_as.extension().and_then(|x| x.to_str()) {
            Some("yaml") => serde_yaml::to_writer(&f, self).expect("error writing checkpoint?!"),
            Some("json") => serde_json::to_writer(&f, self).expect("error writing checkpoint?!"),
            Some("cbor") => serde_cbor::to_writer(&f, self).expect("error writing checkpoint?!"),
            _ => panic!("I don't know how to create file {:?}", self.save_as),
        }
    }

    fn e_to_idx(&self, energy: Energy) -> usize {
        bisect(&self.energy_boundaries, energy)
    }

    /// Run a simulation
    pub fn run_once(&mut self) {
        self.moves += 1;
        self.rejected_moves += 1;
        let e1 = self.system.energy();
        let i1 = self.e_to_idx(e1);
        self.histogram[i1] += 1;
        self.total_energy[i1] += e1;
        let distance = self.system.max_size() * self.rng.gen_range(-10.0, 0.0f64).exp();
        if let Some(e2) = self.system.plan_move(&mut self.rng, distance) {
            let i2 = self.e_to_idx(e2);
            let lnw1 = self.lnw[i1];
            let lnw2 = self.lnw[i2];
            if lnw2 <= lnw1 || self.rng.gen::<f64>() < (lnw1 - lnw2).exp() {
                self.system.confirm();
                self.rejected_moves -= 1;
            }
        }

        let movie_time = self.movie.shall_i_save(self.moves);
        if movie_time {
            self.movie.save_frame(&self.save_as, self.moves, &self);
        }
        if self.report.am_all_done(self.moves, self.moves) {
            self.checkpoint();
            println!("All done!");
            ::std::process::exit(0);
        }
        if self.save.shall_i_save(self.moves) || movie_time {
            self.checkpoint();
        }
    }
}

fn bisect(b: &[Energy], v: Energy) -> usize {
    match b.binary_search_by(|&x| {
        if v < x {
            std::cmp::Ordering::Less
        } else {
            std::cmp::Ordering::Greater
        }
    }) {
        Err(i) => i,
        Ok(_) => unreachable!(),
    }
}

#[test]
fn bisect_ok() {
    let b = &[
        Energy::new(5.0),
        Energy::new(3.0),
        Energy::new(1.0),
        Energy::new(-1.0),
        Energy::new(-2.0),
    ];
    assert_eq!(bisect(b, Energy::new(6.0)), 0);
    assert_eq!(bisect(b, Energy::new(4.0)), 1);
    assert_eq!(bisect(b, Energy::new(2.0)), 2);
    assert_eq!(bisect(b, Energy::new(0.0)), 3);
    assert_eq!(bisect(b, Energy::new(-1.1)), 4);
    assert_eq!(bisect(b, Energy::new(-2.2)), 5);
    let mut b = Vec::new();
    for i in 0..138 {
        b.push(Energy::new(-i as f64));
    }
    for i in 0..139 {
        println!("looking for {}", i);
        assert_eq!(bisect(&b, Energy::new(-(i as f64) + 0.5)), i);
    }
}
