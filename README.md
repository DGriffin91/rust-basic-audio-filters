# Audio Filters
A collection of filters for real-time audio processing

### Filter Types

- [x] Bell
- [x] Low Pass
- [x] High Pass
- [x] Low Shelf
- [x] High Shelf
- [x] Notch
- [x] Band Pass
- [x] All Pass

### Filter Features

- [x] Bode Plot (phase & amplitude)
- [x] 1st and 2nd order filter primitives
- [x] Virtual analog (VA) State Variable Filters (SVF) for both 1st & 2nd order IIR.
- [x] Minimum Phase IIR Mode

```rust
let fs = 48000.0;
let f0 = 1000.0;
let db_gain = 6.0;
let q_value = 1.0;

let coeffs = IIR2Coefficients::highshelf(f0, db_gain, q_value, fs);

let mut filter = IIR2::from(coeffs);

for i in 0..1000 {
    audio[i] = filter.process(audio[i]);
}
```

