use std::f32::consts::{PI, TAU};

use num_complex::Complex;

pub fn get_z(f_hz: f32, fs: f32) -> Complex<f32> {
    let z = -TAU * f_hz / fs;
    z.cos() + z.sin() * Complex::<f32>::new(0.0, 1.0)
}

#[derive(Copy, Clone, Debug)]
pub struct IIR1Coefficients {
    pub a: f32,
    pub g: f32,
    pub a1: f32,
    pub m0: f32,
    pub m1: f32,
    pub fs: f32,
}

impl IIR1Coefficients {
    pub fn get_bode_sample(self, f_hz: f32, fs: f32) -> Complex<f32> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.

        let z = -TAU * f_hz / fs;
        let z = z.cos() + z.sin() * Complex::<f32>::new(0.0, 1.0);

        let denominator = self.g + z * (self.g - 1.0) + 1.0;

        let y = self.m0 + (self.m1 * self.g * (z + 1.0)) / denominator;

        y
    }

    pub fn empty() -> IIR1Coefficients {
        IIR1Coefficients {
            a: 0.0,
            g: 0.0,
            a1: 0.0,
            m0: 0.0,
            m1: 0.0,
            fs: 0.0,
        }
    }

    pub fn lowpass(f0: f32, _db_gain: f32, fs: f32) -> IIR1Coefficients {
        let f0 = f0.min(fs * 0.0_5);
        let a = 1.0;
        let g = (PI * f0 / fs).tan();
        let a1 = g / (1.0 + g);
        let m0 = 0.0;
        let m1 = 1.0;
        IIR1Coefficients {
            a,
            g,
            a1,
            m0,
            m1,
            fs,
        }
    }

    pub fn highpass(f0: f32, _db_gain: f32, fs: f32) -> IIR1Coefficients {
        let f0 = f0.min(fs * 0.0_5);
        let a = 1.0;
        let g = (PI * f0 / fs).tan();
        let a1 = g / (1.0 + g);
        let m0 = 1.0;
        let m1 = -1.0;
        IIR1Coefficients {
            a,
            g,
            a1,
            m0,
            m1,
            fs,
        }
    }

    pub fn allpass(f0: f32, _db_gain: f32, fs: f32) -> IIR1Coefficients {
        let f0 = f0.min(fs * 0.0_5);
        let a = 1.0;
        let g = (PI * f0 / fs).tan();
        let a1 = g / (1.0 + g);
        let m0 = 1.0;
        let m1 = -2.0;
        IIR1Coefficients {
            a,
            g,
            a1,
            m0,
            m1,
            fs,
        }
    }

    pub fn lowshelf(f0: f32, db_gain: f32, fs: f32) -> IIR1Coefficients {
        let f0 = f0.min(fs * 0.0_5);
        let a = 1.0f32.powf(db_gain / 20.0);
        let g = (PI * f0 / fs).tan() / (a).sqrt();
        let a1 = g / (1.0 + g);
        let m0 = 1.0;
        let m1 = a - 1.0;
        IIR1Coefficients {
            a,
            g,
            a1,
            m0,
            m1,
            fs,
        }
    }

    pub fn highshelf(f0: f32, db_gain: f32, fs: f32) -> IIR1Coefficients {
        let f0 = f0.min(fs * 0.0_5);
        let a = 1.0f32.powf(db_gain / 20.0);
        let g = (PI * f0 / fs).tan() * (a).sqrt();
        let a1 = g / (1.0 + g);
        let m0 = a;
        let m1 = 1.0 - a;
        IIR1Coefficients {
            a,
            g,
            a1,
            m0,
            m1,
            fs,
        }
    }
}

/// Internal states and coefficients of the SVF form
#[derive(Copy, Clone, Debug)]
pub struct IIR1 {
    ic1eq: f32,
    pub coeffs: IIR1Coefficients,
}

impl IIR1 {
    /// Creates a SVF from a set of filter coefficients
    pub fn from(coefficients: IIR1Coefficients) -> Self {
        IIR1 {
            ic1eq: 0.0,
            coeffs: coefficients,
        }
    }

    pub fn process(&mut self, input: f32) -> f32 {
        let v1 = self.coeffs.a1 * (input - self.ic1eq);
        let v2 = v1 + self.ic1eq;
        self.ic1eq = v2 + v1;

        self.coeffs.m0 * input + self.coeffs.m1 * v2
    }

    pub fn update_coefficients(&mut self, new_coefficients: IIR1Coefficients) {
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
    fn test_iir1() {
        let mut audio: Vec<f32> = (0..1000).map(|x| rand(x as f32)).collect();

        let fs = 48000.0;
        let f0 = 1000.0;
        let db_gain = 6.0;

        let coeffs = IIR1Coefficients::highshelf(f0, db_gain, fs);

        let mut filter = IIR1::from(coeffs);

        for i in 0..1000 {
            audio[i] = filter.process(audio[i]);
        }

        assert_eq!(audio[500], -0.12890625)
    }
}
