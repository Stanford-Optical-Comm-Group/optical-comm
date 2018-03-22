# Processing

This folder contains functions for signal processing and impairment mitigation.

## Contents 

Below is a list of classes and functions in this folder. Some examples are shown, but more detailed documentation can be found in the files.

### Equalization

- [equalize.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/processing/equalize.m): **Design and implementation of linear equalizer**

Designs and implement a linear equalizer. Design can be done adaptively or not.

Function call to design linear equalizer for signal `yt`.
```
[yd, eq] = equalize(eq, yt, Hch, ModFormat, sim, verbose)
```

- [fse_cma.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/processing/fse_cma.m): **Fractionally spaced equalization (FSE) using constant modulus algorithm (CMA)**

Fractionally spaced equalization (FSE) using constant modulus algorithm (CMA). Equalization is done in time domain.

This function can use either four or two filters in the 2 x 2 MIMO equalizer. The two-filter option assumes negligible DGD and is useful to halve the number of operations in DSP.

Function call:
```
[Y, W, MSE] = fse_cma(X, eq, verbose)
```

- [lms_td_fse.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/processing/lms_td_fse.m): **Adaptive time-domain fractionally spaced linear equalizer adapted using LMS algorithm**

This function assumes that there is no phase noise.

This function can use either four or two filters in the 2 x 2 MIMO equalizer. The two-filter option assumes negligible DGD and is useful to halve the number of operations in DSP.

Function call:
```
[Y, W, MSE] = lms_td_fse(X, eq, verbose)
```

### Carrier phase recovery

- [feedforward_cpr.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/processing/feedforward_cpr.m): **Feedforward carrier phase recovery**

Simulates carrier recovery for QPSK in two polarizations.

Function call:
```
[xhat, thetahat, Delta] = feedforward_cpr(y, CPR, ModFormat, verbose)
```

- [phase_estimation.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/processing/phase_estimation.m): **Phase estimation**

Adaptively calculates and correct for phase offset.

Function call:
```
[Y, phi] = phase_estimation(X, eq, M, verbose) 
```

- [qpsk_freq_rec.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/processing/qpsk_freq_rec.m): **QPSK frequency recovery**

Frequency recovery algorithm for QPSK.

Function call:
```
[Y, Df] = qpsk_freq_rec(X, FreqRec,  verbose)
```