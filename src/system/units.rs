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

impl<A> UNITS<f64, A> {
    /// Format the number in a nice way
    pub fn pretty(&self) -> crate::prettyfloat::PrettyFloat {
        crate::prettyfloat::PrettyFloat(*self.value_unsafe())
    }
}

pub use self::f64consts::*;

impl<A> std::iter::Sum<Self> for UNITS<f64, A> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        iter.fold(UNITS::<f64,A>::new(0.0), |a, b| a+b)
    }
}