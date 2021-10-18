use std::f64::consts::{PI, TAU};

use num_complex::Complex;

#[derive(Copy, Clone, Debug)]
pub struct IIR1Coefficients {
    pub a: f64,
    pub g: f64,
    pub a1: f64,
    pub m0: f64,
    pub m1: f64,
}

impl IIR1Coefficients {
    #[inline]
    pub fn get_bode_sample(self, frequency_hz: f64, sample_rate_hz: f64) -> Complex<f64> {
        //Use y.norm() for amplitude and y.arg().to_degrees() for phase. Add to combine phase.

        let z = -TAU * frequency_hz / sample_rate_hz;
        let z = z.cos() + z.sin() * Complex::<f64>::new(0.0, 1.0);

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
        }
    }

    #[inline]
    pub fn lowpass(cutoff_hz: f64, _gain_db: f64, sample_rate_hz: f64) -> IIR1Coefficients {
        let cutoff_hz = cutoff_hz.min(sample_rate_hz * 0.5);
        let a = 1.0;
        let g = (PI * cutoff_hz / sample_rate_hz).tan();
        let a1 = g / (1.0 + g);
        let m0 = 0.0;
        let m1 = 1.0;
        IIR1Coefficients { a, g, a1, m0, m1 }
    }

    #[inline]
    pub fn highpass(cutoff_hz: f64, _gain_db: f64, sample_rate_hz: f64) -> IIR1Coefficients {
        let cutoff_hz = cutoff_hz.min(sample_rate_hz * 0.5);
        let a = 1.0;
        let g = (PI * cutoff_hz / sample_rate_hz).tan();
        let a1 = g / (1.0 + g);
        let m0 = 1.0;
        let m1 = -1.0;
        IIR1Coefficients { a, g, a1, m0, m1 }
    }

    #[inline]
    pub fn allpass(cutoff_hz: f64, _gain_db: f64, sample_rate_hz: f64) -> IIR1Coefficients {
        let cutoff_hz = cutoff_hz.min(sample_rate_hz * 0.5);
        let a = 1.0;
        let g = (PI * cutoff_hz / sample_rate_hz).tan();
        let a1 = g / (1.0 + g);
        let m0 = 1.0;
        let m1 = -2.0;
        IIR1Coefficients { a, g, a1, m0, m1 }
    }

    #[inline]
    pub fn lowshelf(cutoff_hz: f64, gain_db: f64, sample_rate_hz: f64) -> IIR1Coefficients {
        let cutoff_hz = cutoff_hz.min(sample_rate_hz * 0.5);
        let a = 10.0f64.powf(gain_db / 20.0);
        let g = (PI * cutoff_hz / sample_rate_hz).tan() / (a).sqrt();
        let a1 = g / (1.0 + g);
        let m0 = 1.0;
        let m1 = a - 1.0;
        IIR1Coefficients { a, g, a1, m0, m1 }
    }

    #[inline]
    pub fn highshelf(cutoff_hz: f64, gain_db: f64, sample_rate_hz: f64) -> IIR1Coefficients {
        let cutoff_hz = cutoff_hz.min(sample_rate_hz * 0.5);
        let a = 10.0f64.powf(gain_db / 20.0);
        let g = (PI * cutoff_hz / sample_rate_hz).tan() * (a).sqrt();
        let a1 = g / (1.0 + g);
        let m0 = a;
        let m1 = 1.0 - a;
        IIR1Coefficients { a, g, a1, m0, m1 }
    }
}

/// Internal states and coefficients of the SVF form
#[derive(Copy, Clone, Debug)]
pub struct IIR1 {
    ic1eq: f64,
    pub coeffs: IIR1Coefficients,
}

impl IIR1 {
    /// Creates a SVF from a set of filter coefficients
    #[inline]
    pub fn from(coefficients: IIR1Coefficients) -> Self {
        IIR1 {
            ic1eq: 0.0,
            coeffs: coefficients,
        }
    }

    #[inline]
    pub fn process(&mut self, input_sample: f64) -> f64 {
        let v1 = self.coeffs.a1 * (input_sample - self.ic1eq);
        let v2 = v1 + self.ic1eq;
        self.ic1eq = v2 + v1;

        self.coeffs.m0 * input_sample + self.coeffs.m1 * v2
    }

    #[inline]
    pub fn update(&mut self, new_coefficients: IIR1Coefficients) {
        self.coeffs = new_coefficients;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rand(x: f64) -> f64 {
        ((x * 12.9898).sin() * 43758.5453).fract()
    }

    #[test]
    fn test_iir1() {
        let mut audio: Vec<f64> = (0..1000).map(|x| rand(x as f64)).collect();

        let sample_rate_hz = 48000.0;
        let cutoff_hz = 1000.0;
        let gain_db = 6.0;

        let coeffs = IIR1Coefficients::highshelf(cutoff_hz, gain_db, sample_rate_hz);

        let mut filter = IIR1::from(coeffs);

        for i in 0..1000 {
            audio[i] = filter.process(audio[i]);
        }

        assert_eq!(audio[500], -0.9407069884176158)
    }
}
