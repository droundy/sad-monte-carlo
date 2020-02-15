//! Energy/distance units

make_units! {
    UNITS;
    ONE: Unitless;

    base {
        EPSILON: Energy, "e", Energy;
        SIGMA: Distance, "d", Length;
    }

    derived {
        AREA: Area = (Distance * Distance), Area;
        VOLUME: Volume = (Distance * Distance * Distance), Volume;
        DENSITY: Density = (Unitless / Distance / Distance / Distance);
        PRESSURE: Pressure = (Energy / Distance / Distance / Distance);
        FORCE: Force = (Energy / Distance);
        ENERGY_SQUARED: EnergySquared = (Energy*Energy);
        PER_ENERGY: PerEnergy = (Unitless / Energy);
    }

    constants {
        PI: Unitless = consts::PI;
        R: Distance = 0.5;
    }

    fmt = true;
}

impl_serde!(UNITS);
impl_auto_args!(UNITS);

impl Unitless<f64> {
    /// Format the number in a nice way
    pub fn pretty(&self) -> crate::prettyfloat::PrettyFloat {
        crate::prettyfloat::PrettyFloat(*self.value())
    }
}
impl Energy<f64> {
    /// Format the number in a nice way
    pub fn pretty(&self) -> crate::prettyfloat::PrettyFloat {
        crate::prettyfloat::PrettyFloat(*(*self/EPSILON).value())
    }
}

pub use self::f64consts::*;
