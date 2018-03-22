# Measurement

This folder contains functions and classes use for making measurements. 

Files with capitalized letters are classes, while files in lower case are functions.

## Contents 

Below is a list of classes and functions in this folder. Some examples are shown, but more detailed documentation can be found in the files.

### Electronics

- [eyediagram.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/measurement/eyediagram.m): **Eye diagram**

Simple eye-diagram plot. Matlab communications toolbox offers more sophisticated options

Function call to plot eye-diagram of signal `x` whose symbol period is `n`
```
eyediagram(x, n)
```

### Optics
- [OSA.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/measurement/OSA.m): **Optical Spectrum Analyzer**

Class for optical spectrum analyzer. This class produces a plot of the measured optical spectrum and also estimaters OSNR

Constructor for an OSA with resolution 0.1 nm
```
Osa = OSA(0.1e-9);
```

Show optical spectrum of signal `E` at 1550e-9
```
Osa.spectrum(E, 1550e-9, frequency_vector, 'true');
```

Estimate OSNR optical spectrum of signal `E` at 1550e-9
```
Osa = OSA(0.4e-9);
[OSNRdB, OSNRdB01nm] = Osa.estimate_osnr(E, 1550e-9, frequency_vector, 'true')
```
It's usually necessary to use an OSA with resolution higher than 0.1 nm in order to measure all the signal power (e.g., 0.4 nm). The function `estimate_osnr` automatically converts the measurement to the 0.1 nm reference bandwidth (`OSNRdB01nm`).

- [power_meter.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/measurement/power_meter.m): **Power meter**

Simply measures optical power. Signal can be in multiple dimensions (modes or polarizations)

Function call to measure optical power in dBm and Watts of signal `x`
```
[PdBm, PW] = power_meter(x)
```