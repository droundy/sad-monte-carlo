//! Energy/distance units

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
    }

    fmt = true;
}

impl_serde!(UNITS);

pub use self::f64consts::*;
