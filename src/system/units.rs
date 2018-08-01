//! Energy/distance units

use clapme::ClapMe;

make_units! {
    UNITS;
    ONE: Unitless;

    base {
        EPSILON: Energy, "epsilon", Energy;
        SIGMA: Distance, "sigma", Length;
    }

    derived {
        AREA: Area = (Distance * Distance), Area;
        VOLUME: Volume = (Distance * Distance * Distance), Volume;
        DENSITY: Density = (Unitless / Distance / Distance / Distance);
        PRESSURE: Pressure = (Energy / Distance / Distance / Distance);
        FORCE: Force = (Energy / Distance);
    }

    constants {
        PI: Unitless = consts::PI;
        R: Distance = 0.5;
    }

    fmt = true;
}

impl_serde!(UNITS);

impl<V: ClapMe,U> ClapMe for UNITS<V,U> {
    fn with_clap<T>(info: ::clapme::ArgInfo, app: ::clapme::clap::App,
                    f: impl FnOnce(::clapme::clap::App) -> T) -> T {
        V::with_clap(info, app, f)
    }
    fn from_clap(name: &str, matches: &::clapme::clap::ArgMatches) -> Option<Self> {
        V::from_clap(name, matches).map(|v| UNITS {
            value_unsafe: v,
            _marker: marker::PhantomData,
        })
    }
    fn requires_flags(name: &str) -> Vec<String> {
        V::requires_flags(name)
    }
}

pub use self::f64consts::*;
