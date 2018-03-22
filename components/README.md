# Components

This folder contains components and devices. 

Files with capitalized letters are classes, while files in lower case are functions.

## Contents 

Below is a list of classes and functions in this folder. Some examples are shown, but more detailed documentation can be found in the files.

### Electronics

- [adc.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/components/adc.m): **Analog-to-digital converter**
Filter, downsample, and quantize signal `x`
Function call for converting analog signal `x` to digitized signal `ys` using ADC parameters `ADC` and simulation parameters `sim`. See code for details
```
[ys, varQ, xa] = adc(x, ADC, sim)
```

- [dac.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/components/dac.m): **Digital-to-analog converter**
Quantize, sample and hold with ZOH, and analog filtering
Function call for converting digital signal `x` to analog signal `xa` using DAC parameters `DAC` and simulation parameters `sim`. See code for details
```
[xa, xqzoh, xq] = dac(x, DAC, sim, verbose)
```

- [Filter.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/components/Filter.m): **Filter**
This class can be used as a filter within a loop. It can take one value at each iteration. This is useful for simulating a filter within a loop. When used this way, the filter is implemented as a direct form I. If a vector is passed as a parameter, then the filtering operation will be performed using Matlab's `filter` function (transposed direct form I). If the flag `removeGroupDelay` is `true`, the filter will be implemented in the frequency domain.

Constructor options:
..* filter specifications are passed: type, order, normalized cutoff frequency, fs (optional)
..* filter struct designed with design_filter.m
..* filter is already known: passing num, den, fs (optional)

Example of sixth-order Butterworth filter with cutoff frequency 0.5
```
F = Filter('butter', 6, 0.5);

filt = design_filter('butter', 6, 0.5);
F = Filter(filt)

[num, den] = butter(6, 0.5);
F = Filter(num, den);
```

### Photodiodes
- [APD.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/components/APD.m): **Avalanche photodiode**
Class for modeling avalanche photodiode.

Constructor:
```
Apd = APD(GaindB, ka, BW, R, Id)
```

Example of APD with gain 10 dB, impact ionization factor (`ka`) 0.1, low-gain bandwidth (`BW(1)`) 30 GHz, gain-bandwidth product (`BW(2)`) 300 GHz, responsitivity (`R`) 0.7, and dark current (`Id`) 10 nA

```
Apd = APD(10, 0.1, [30 300], 0.7, 10e-9);
current = Apd.detect(electric_field, sampling_frequency, 'Gaussian', N0);
```
In this case the APD method `detect` also adds thermal noise with one-sided PSD `N0`.

- [PIN.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/components/PIN.m): **Positive-intrinsic-negative photodiode**
This class is a daughter class of class [APD.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/components/APD.m).

Constructor:
```
PD = PIN(R, Id, BW)
```

Example of PIN responsitivity (`R`) 1, bandwidth `BW` 30 GHz, and dark current (`Id`) 10 nA.

```
PD = PIN(1, 30, 10e-9);
current = PD.detect(electric_field, sampling_frequency, 'Gaussian', N0);
```
In this case the PIN method `detect` also adds thermal noise with one-sided PSD `N0`.


### Fibers and waveguides

- [EDF.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/components/EDF.m): **Single-mode erbium-doped fiber**
Class for modeling of Erbium-doped fiber physics using analytical, semi-analytical, and numerical methods.

Constructor:
```
E = EDF(L, type);
```

- [Fiber.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/components/EDF.m): **Single-mode fiber**
Class for modeling single-mode fiber

Constructor:
```
SMF = Fiber(L, att, D)
```
`L` is the fiber length in m, `att` is a function handle of the fiber attenuation in dB/km as a function of the wavelength (default `att = @(l) 0`, and `D` is a function handle of fiber chromatic dispersion (s/m^2) as a function of the wavelength (default SMF28 dispersion, see code for details).

Example of linear propagation over 100 km of SMF28 fiber for a signal at 1550 nm.
```
SMF28 = Fiber(100e3, @(l) 0.2);
Eout = SMF28.linear_propagation(Ein, frequency_vector, 1550e-9);
```

### Laser, modulators, and amplifiers

- [Laser.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/components/Laser.m): **Single-mode Laser**
Modeling includes intensity noise, phase noise, chirp, frequency offset, and frequency response for directly modulated lasers

Constructor options:
```
LO = Laser(lambda, PdBm, RIN, linewidth, freqOffset);
```

- [mzm.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/components/mzm.m): **Mach-Zehnder modulator*

Mach-Zehnder modulator including bandwidth limitation, nonlinear transfer function, and imperfect biasing and clipping.

Function call for driving signal `Vin`, input electric field `Ein`, and modulator parameters `Mod`. Among other things, `Mod` defines the MZM frequency response. If the MZM is traveling wave device its frequency response is limited by loss and velocity mismatch, whereas if the MZM is a lumped device its frequency response is characterized by parasitics. 
```
Eout = mzm(Ein, Vin, Mod)
```

Type of MZ modulators:
..* __Dual-pol I/Q modulator__ if `Vin` is 2 x N and complex
..* __Single-pol I/Q modulator__ if `Vin` is 1 x N and complex
..* __Intensity modulator__ if `Vin` is 1 x N and real

- [OpticalAmplifier.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/components/OpticalAmplifier.m): **Single-mode Optical Amplifier**

Class for modeling optical amplifier. Amplifier is modeled based on gain, noise figure. Hence, to first order, this class can be used to model both EDFAs and SOAs.

Constructor:
```
Amp = OpticalAmplifier(Operation, param, NoiseFigure, Wavelength);
```

Two operations are possible: `ConstantOutputPower` and `ConstantGain`. If `Operation = ConstantOutputPower`, `param` is interpreted as the desired output power. Similarly, if `Operation = ConstantGain`, `param` is interpreted as the desired gain.


[Eout1, Eout2] = power_splitter(Ein1, Ein2, splittingRatio, phaseShiftDeg, timeDelay, f)

### Misc
- [attenuator.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/components/attenuator.m): **Attenuator**
Attenuate signal `x` by `attdB` dB

```
x = attenuator(x, attdB);
```

- [power_splitter.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/components/attenuator.m): **Power splitter**

Function call:
```
[Eout1, Eout2] = power_splitter(Ein1, Ein2, splittingRatio, phaseShiftDeg, timeDelay, f)
```

Example of ideal 3-dB power splitter:

```
[Eout1, Eout2] = power_splitter(Ein1, Ein2, 0.5, 0);
```

This results in the transfer matrix: `[Eout1; Eout2] = 1/sqrt(2)*[1 j; j 1][Ein1; Ein2]`
