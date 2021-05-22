use std::f32::consts::{PI, TAU};

use num_complex::Complex;

#[derive(Copy, Clone, Debug)]
pub struct IIR2Coefficients {
    pub a: f32,
    pub g: f32,
    pub gpow2: f32,
    pub k: f32,
    pub a1: f32,
    pub a2: f32,
    pub a3: f32,
    pub m0: f32,
    pub m1: f32,
    pub m2: f32,
    pub fs: f32,
}

impl IIR2Coefficients {
    pub fn get_bode_sample(self, f_hz: f32, fs: f32) -> Complex<f32> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.

        let z = -TAU * f_hz / fs;
        let z = z.cos() + z.sin() * Complex::<f32>::new(0.0, 1.0);
        let zpow2 = z * z;

        let denominator = (self.gpow2 + self.g * self.k + 1.0)
            + 2.0 * (self.gpow2 - 1.0) * z
            + (self.gpow2 - self.g * self.k + 1.0) * zpow2;

        let y = self.m0
            + (self.m1 * self.g * (1.0 - zpow2) + self.m2 * self.gpow2 * (1.0 + 2.0 * z + zpow2))
                / denominator;

        y
    }

    pub fn lowpass(f0: f32, q_value: f32, _db_gain: f32, fs: f32) -> IIR2Coefficients {
        let f0 = f0.min(fs * 0.5);
        let a = 1.0;
        let g = (PI * f0 / fs).tan();
        let k = 1.0 / q_value;
        let a1 = 1.0 / (1.0 + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;
        let m0 = 0.0;
        let m1 = 0.0;
        let m2 = 1.0;
        IIR2Coefficients {
            a,
            g,
            gpow2: g * g,
            k,
            a1,
            a2,
            a3,
            m0,
            m1,
            m2,
            fs,
        }
    }
    pub fn highpass(f0: f32, q_value: f32, _db_gain: f32, fs: f32) -> IIR2Coefficients {
        let f0 = f0.min(fs * 0.5);
        let a = 1.0;
        let g = (PI * f0 / fs).tan();
        let k = 1.0 / q_value;
        let a1 = 1.0 / (1.0 + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;
        let m0 = 1.0;
        let m1 = -k;
        let m2 = -1.0;
        IIR2Coefficients {
            a,
            g,
            gpow2: g * g,
            k,
            a1,
            a2,
            a3,
            m0,
            m1,
            m2,
            fs,
        }
    }
    pub fn bandpass(f0: f32, q_value: f32, _db_gain: f32, fs: f32) -> IIR2Coefficients {
        let f0 = f0.min(fs * 0.5);
        let a = 1.0;
        let g = (PI * f0 / fs).tan();
        let k = 1.0 / q_value;
        let a1 = 1.0 / (1.0 + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;
        let m0 = 0.0;
        let m1 = 1.0;
        let m2 = 0.0;
        IIR2Coefficients {
            a,
            g,
            gpow2: g * g,
            k,
            a1,
            a2,
            a3,
            m0,
            m1,
            m2,
            fs,
        }
    }
    pub fn notch(f0: f32, q_value: f32, _db_gain: f32, fs: f32) -> IIR2Coefficients {
        let f0 = f0.min(fs * 0.5);
        let a = 1.0;
        let g = (PI * f0 / fs).tan();
        let k = 1.0 / q_value;
        let a1 = 1.0 / (1.0 + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;
        let m0 = 1.0;
        let m1 = -k;
        let m2 = 0.0;
        IIR2Coefficients {
            a,
            g,
            gpow2: g * g,
            k,
            a1,
            a2,
            a3,
            m0,
            m1,
            m2,
            fs,
        }
    }
    pub fn allpass(f0: f32, q_value: f32, _db_gain: f32, fs: f32) -> IIR2Coefficients {
        let f0 = f0.min(fs * 0.5);
        let a = 1.0;
        let g = (PI * f0 / fs).tan();
        let k = 1.0 / q_value;
        let a1 = 1.0 / (1.0 + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;
        let m0 = 1.0;
        let m1 = -2.0 * k;
        let m2 = 0.0;
        IIR2Coefficients {
            a,
            g,
            gpow2: g * g,
            k,
            a1,
            a2,
            a3,
            m0,
            m1,
            m2,
            fs,
        }
    }
    pub fn lowshelf(f0: f32, q_value: f32, db_gain: f32, fs: f32) -> IIR2Coefficients {
        let f0 = f0.min(fs * 0.5);
        let a = 10.0f32.powf(db_gain / 40.0);
        let g = (PI * f0 / fs).tan() / a.sqrt();
        let k = 1.0 / q_value;
        let a1 = 1.0 / (1.0 + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;
        let m0 = 1.0;
        let m1 = k * (a - 1.0);
        let m2 = a * a - 1.0;
        IIR2Coefficients {
            a,
            g,
            gpow2: g * g,
            k,
            a1,
            a2,
            a3,
            m0,
            m1,
            m2,
            fs,
        }
    }
    pub fn highshelf(f0: f32, q_value: f32, db_gain: f32, fs: f32) -> IIR2Coefficients {
        let f0 = f0.min(fs * 0.5);
        let a = 10.0f32.powf(db_gain / 40.0);
        let g = (PI * f0 / fs).tan() * a.sqrt();
        let k = 1.0 / q_value;
        let a1 = 1.0 / (1.0 + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;
        let m0 = a * a;
        let m1 = k * (1.0 - a) * a;
        let m2 = 1.0 - a * a;
        IIR2Coefficients {
            a,
            g,
            gpow2: g * g,
            k,
            a1,
            a2,
            a3,
            m0,
            m1,
            m2,
            fs,
        }
    }
    pub fn bell(f0: f32, q_value: f32, db_gain: f32, fs: f32) -> IIR2Coefficients {
        let f0 = f0.min(fs * 0.5);
        let a = 10.0f32.powf(db_gain / 40.0);
        let g = (PI * f0 / fs).tan();
        let k = 1.0 / (q_value * a);
        let a1 = 1.0 / (1.0 + g * (g + k));
        let a2 = g * a1;
        let a3 = g * a2;
        let m0 = 1.0;
        let m1 = k * (a * a - 1.0);
        let m2 = 0.0;
        IIR2Coefficients {
            a,
            g,
            gpow2: g * g,
            k,
            a1,
            a2,
            a3,
            m0,
            m1,
            m2,
            fs,
        }
    }
}

/// Internal states and coefficients of the SVF form
#[derive(Copy, Clone, Debug)]
pub struct IIR2 {
    ic1eq: f32,
    ic2eq: f32,
    pub coeffs: IIR2Coefficients,
}

impl IIR2 {
    /// Creates a SVF from a set of filter coefficients
    pub fn from(coefficients: IIR2Coefficients) -> Self {
        IIR2 {
            ic1eq: 0.0,
            ic2eq: 0.0,
            coeffs: coefficients,
        }
    }

    pub fn process(&mut self, input: f32) -> f32 {
        let v3 = input - self.ic2eq;
        let v1 = self.coeffs.a1 * self.ic1eq + self.coeffs.a2 * v3;
        let v2 = self.ic2eq + self.coeffs.a2 * self.ic1eq + self.coeffs.a3 * v3;
        self.ic1eq = 2.0 * v1 - self.ic1eq;
        self.ic2eq = 2.0 * v2 - self.ic2eq;

        self.coeffs.m0 * input + self.coeffs.m1 * v1 + self.coeffs.m2 * v2
    }

    pub fn update_coefficients(&mut self, new_coefficients: IIR2Coefficients) {
        self.coeffs = new_coefficients;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rand(x: f32) -> f32 {
        ((x * 12.9898).sin() * 43758.5453).fract()
    }

    #[test]
    fn test_iir2() {
        let mut audio: Vec<f32> = (0..1000).map(|x| rand(x as f32)).collect();

        let fs = 48000.0;
        let f0 = 1000.0;
        let db_gain = 6.0;
        let q_value = 1.0;

        let coeffs = IIR2Coefficients::highshelf(f0, q_value, db_gain, fs);

        let mut filter = IIR2::from(coeffs);

        for i in 0..1000 {
            audio[i] = filter.process(audio[i]);
        }

        assert_eq!(audio[500], -0.5090322)
    }
}
