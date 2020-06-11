//! This module attempts to provide a %g formatter.

use std::fmt::{Alignment, Display, Formatter, Result};

/// Wrap this type around an `f64` in order to print it nicely.
pub struct PrettyFloat(pub f64);

fn n_decimals(value: f64, digits: usize) -> usize {
    if value == 0. {
        return 0;
    }
    let log10 = value.abs().log10();
    if log10 > digits as f64 {
        0
    } else {
        (digits as f64 - log10) as usize + 1
    }
}
#[test]
fn test_n_decimals() {
    assert_eq!(n_decimals(10.1, 1), 0);
    assert_eq!(n_decimals(10.1, 2), 1);
    assert_eq!(n_decimals(0.02, 1), 3);
    assert_eq!(n_decimals(0.002, 1), 4);
}

fn cut_trailing(mut x: String) -> String {
    if x.contains(".") {
        if x.contains("e") {
            let parts: Vec<_> = x.splitn(2, 'e').collect();
            let mut x = parts[0].to_string();
            while x.ends_with('0') {
                x.pop();
            }
            if x.ends_with('.') {
                x.pop();
            }
            x.push('e');
            x.push_str(parts[1]);
            x
        } else {
            while x.ends_with('0') {
                x.pop();
            }
            if x.ends_with('.') {
                x.pop();
            }
            x
        }
    } else {
        x
    }
}

impl Display for PrettyFloat {
    fn fmt(&self, f: &mut Formatter) -> Result {
        if let Some(precision) = f.precision() {
            let ndec = n_decimals(self.0, precision);
            let options = &[
                cut_trailing(format!("{}", self.0)),
                cut_trailing(format!("{:e}", self.0)),
                cut_trailing(format!("{:.*}", ndec, self.0)),
                cut_trailing(format!("{:.*e}", precision, self.0)),
            ];
            let s = options.into_iter().min_by_key(|s| s.len()).unwrap();
            if let Some(width) = f.width() {
                // If we received a width, we use it
                match f.align() {
                    Some(Alignment::Left) => write!(f, "{:<width$}", s, width = width),
                    Some(Alignment::Right) => write!(f, "{:>width$}", s, width = width),
                    Some(Alignment::Center) => write!(f, "{:^width$}", s, width = width),
                    None => write!(f, "{:width$}", s, width = width),
                }
            } else {
                f.write_str(s)
            }
        } else {
            let options = &[format!("{}", self.0), format!("{:e}", self.0)];
            let s = options.into_iter().min_by_key(|s| s.len()).unwrap();
            if let Some(width) = f.width() {
                // If we received a width, we use it
                match f.align() {
                    Some(Alignment::Left) => write!(f, "{:<width$}", s, width = width),
                    Some(Alignment::Right) => write!(f, "{:>width$}", s, width = width),
                    Some(Alignment::Center) => write!(f, "{:^width$}", s, width = width),
                    None => write!(f, "{:width$}", s, width = width),
                }
            } else {
                f.write_str(s)
            }
        }
    }
}

#[test]
fn known_numbers() {
    println!("0.1 default");
    assert_eq!(&format!("{}", PrettyFloat(0.1)), "0.1");
    println!("0.1 prec=1");
    assert_eq!(&format!("{:.1}", PrettyFloat(0.1)), "0.1");
    println!("0.1 prec=2");
    assert_eq!(&format!("{:.2}", PrettyFloat(0.1)), "0.1");
    println!("0.1 prec=3");
    assert_eq!(&format!("{:.3}", PrettyFloat(0.1)), "0.1");
    println!("0.1 prec=4");
    assert_eq!(&format!("{:.4}", PrettyFloat(0.1)), "0.1");

    println!("0.1 prec=1 width width");
    assert_eq!(&format!("{:5}", PrettyFloat(0.1)), "0.1  ");
    assert_eq!(&format!("{:<5}", PrettyFloat(0.1)), "0.1  ");
    assert_eq!(&format!("{:^5}", PrettyFloat(0.1)), " 0.1 ");
    assert_eq!(&format!("{:>5}", PrettyFloat(0.1)), "  0.1");

    println!("10.1");
    assert_eq!(&format!("{}", PrettyFloat(10.1)), "10.1");
    println!("10.1 prec 1");
    assert_eq!(&format!("{:.1}", PrettyFloat(10.1)), "10");
    println!("10.1 prec 2");
    assert_eq!(&format!("{:.2}", PrettyFloat(10.1)), "10.1");
    println!("10.1 prec 3");
    assert_eq!(&format!("{:.3}", PrettyFloat(10.1)), "10.1");
    println!("10.1 prec 4");
    assert_eq!(&format!("{:.4}", PrettyFloat(10.1)), "10.1");

    println!("1.11111111");
    assert_eq!(&format!("{}", PrettyFloat(1.11111111)), "1.11111111");
    println!("1.11111111 prec 1");
    assert_eq!(&format!("{:.1}", PrettyFloat(1.11111111)), "1.1");
    println!("1.11111111 prec 2");
    assert_eq!(&format!("{:.2}", PrettyFloat(1.11111111)), "1.11");
    println!("1.11111111 prec 3");
    assert_eq!(&format!("{:.3}", PrettyFloat(1.11111111)), "1.111");
    println!("1.11111111 prec 4");
    assert_eq!(&format!("{:.4}", PrettyFloat(1.11111111)), "1.1111");

    println!("1.00000000001");
    assert_eq!(&format!("{}", PrettyFloat(1.00000000001)), "1.00000000001");
    println!("1.11111111 prec 1");
    assert_eq!(&format!("{:.1}", PrettyFloat(1.00000000001)), "1");
    println!("1.11111111 prec 2");
    assert_eq!(&format!("{:.2}", PrettyFloat(1.00000000001)), "1");
    println!("1.11111111 prec 3");
    assert_eq!(&format!("{:.3}", PrettyFloat(1.00000000001)), "1");
    println!("1.11111111 prec 4");
    assert_eq!(&format!("{:.4}", PrettyFloat(1.00000000001)), "1");

    println!("1.00000000001e30");
    assert_eq!(
        &format!("{}", PrettyFloat(1.00000000001e30)),
        "1.00000000001e30"
    );
    println!("1.11111111 prec 1");
    assert_eq!(&format!("{:.1}", PrettyFloat(1.00000000001e30)), "1e30");
    println!("1.11111111 prec 2");
    assert_eq!(&format!("{:.2}", PrettyFloat(1.00000000001e30)), "1e30");
    println!("1.11111111 prec 3");
    assert_eq!(&format!("{:.3}", PrettyFloat(1.00000000001e30)), "1e30");
    println!("1.11111111 prec 4");
    assert_eq!(&format!("{:.4}", PrettyFloat(1.00000000001e30)), "1e30");

    println!("1.11111111e200");
    assert_eq!(
        &format!("{}", PrettyFloat(1.11111111e200)),
        "1.11111111e200"
    );
    println!("1.11111111 prec 1");
    assert_eq!(&format!("{:.1}", PrettyFloat(1.11111111e200)), "1.1e200");
    println!("1.11111111 prec 2");
    assert_eq!(&format!("{:.2}", PrettyFloat(1.11111111e200)), "1.11e200");
    println!("1.11111111 prec 3");
    assert_eq!(&format!("{:.3}", PrettyFloat(1.11111111e200)), "1.111e200");
    println!("1.11111111 prec 4");
    assert_eq!(&format!("{:.4}", PrettyFloat(1.11111111e200)), "1.1111e200");

    println!("1e200");
    assert_eq!(&format!("{}", PrettyFloat(1e200)), "1e200");
    println!("1e200 prec 1");
    assert_eq!(&format!("{:.1}", PrettyFloat(1e200)), "1e200");
    println!("1e200 prec 2");
    assert_eq!(&format!("{:.2}", PrettyFloat(1e200)), "1e200");
    println!("1e200 prec 3");
    assert_eq!(&format!("{:.3}", PrettyFloat(1e200)), "1e200");
    println!("1e200 prec 4");
    assert_eq!(&format!("{:.4}", PrettyFloat(1e200)), "1e200");
}

#[test]
fn short_representation_with_prec() {
    for &prec in &[1, 3, 6, 16, 30] {
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
    println!(
        "{:.prec$} {:.prec$} {:.prec$}",
        f,
        f,
        PrettyFloat(f),
        prec = prec
    );
    assert!(format!("{:.*}", prec, PrettyFloat(f)).len() <= format!("{:.*e}", prec, f).len());
    if let Ok(fprec) = f64::from_str(&format!("{:.*}", prec, PrettyFloat(f))) {
        assert!(((fprec - f) / f).abs() < 10_f64.powf(-(prec as f64)));
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
    assert!(format!("{}", PrettyFloat(f)).len() <= format!("{}", f).len());
    assert!(format!("{}", PrettyFloat(f)).len() <= format!("{:e}", f).len());
    assert_eq!(f64::from_str(&format!("{}", PrettyFloat(f))), Ok(f));
}
