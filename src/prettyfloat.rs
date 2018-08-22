//! This module attempts to provide a %g formatter.

use std::fmt::{Display, Formatter, Result};

/// Wrap this type around an `f64` in order to print it nicely.
pub struct PrettyFloat(f64);

fn n_decimals(value: f64, digits: usize) -> usize {
    let log10 = value.abs().log10();
    if log10 > digits as f64 {
        0
    } else {
        (digits as f64 - log10) as usize + 1
    }
}

impl Display for PrettyFloat {
    fn fmt(&self, f: &mut Formatter) -> Result {
        if let Some(precision) = f.precision() {
            let sf = format!("{}", self.0);
            let se = format!("{:e}", self.0);
            let sf_prec = format!("{:.*}", n_decimals(self.0, precision), self.0);
            let se_prec = format!("{:.*}", precision, self.0);
            if sf.len() < se.len() {
                if sf.len() < sf_prec.len() {
                    if sf.len() < se_prec.len() {
                        write!(f, "{}", sf)
                    } else {
                        write!(f, "{:e}", self.0)
                    }
                } else {
                    if sf_prec.len() < se_prec.len() {
                        write!(f, "{}", self.0)
                    } else {
                        write!(f, "{:.*e}", precision, self.0)
                    }
                }
            } else {
                if se.len() < sf_prec.len() {
                    if se.len() < se_prec.len() {
                        write!(f, "{}", se)
                    } else {
                        write!(f, "{:e}", self.0)
                    }
                } else {
                    if sf_prec.len() < se_prec.len() {
                        write!(f, "{}", self.0)
                    } else {
                        write!(f, "{:.*e}", precision, self.0)
                    }
                }
            }
        } else {
            let sf = format!("{}", self.0);
            let se = format!("{:e}", self.0);
            if se.len() < sf.len() {
                write!(f, "{:e}", self.0)
            } else {
                write!(f, "{}", self.0)
            }
        }
    }
}

#[test]
fn short_representation_with_prec() {
    for &prec in &[1,3,6,16,30] {
        println!("testing with {} digits", prec);
        for &f in &[0.1_f64, 1e-100, 0.1111111111111111] {
            println!("  testing {}", f);
            shortest_with_prec(f, prec);
        }
    }
}

#[cfg(test)]
fn shortest_with_prec(f: f64, prec: usize) {
    use std::str::FromStr;
    println!("{:.prec$} {:.prec$} {:.prec$}", f, f, PrettyFloat(f), prec=prec);
    assert!(format!("{:.*}", prec, PrettyFloat(f)).len() <=
            format!("{:.*e}", prec, f).len());
    if let Ok(fprec) = f64::from_str(&format!("{:.*}", prec, PrettyFloat(f))) {
        assert!(((fprec - f)/f).abs() < 10_f64.powf(-(prec as f64)));
    } else {
        panic!("this is crazy talk!");
    }
}

#[test]
fn short_representation() {
    for &f in &[0.1_f64, 1e-100, 0.1111111111111111] {
        shortest(f);
    }
}

#[cfg(test)]
fn shortest(f: f64) {
    use std::str::FromStr;
    println!("{} {:e} {}", f, f, PrettyFloat(f));
    assert!(format!("{}", PrettyFloat(f)).len() <=
            format!("{}", f).len());
    assert!(format!("{}", PrettyFloat(f)).len() <=
            format!("{:e}", f).len());
    assert_eq!(f64::from_str(&format!("{}", PrettyFloat(f))), Ok(f));
}
