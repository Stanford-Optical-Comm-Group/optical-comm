# Signal

This folder contains classes on modulation formats and signals.

## Contents 

Below is a list of classes and functions in this folder. Some examples are shown, but more detailed documentation can be found in the files.


- [PAM.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/signal/PAM.m): **Pulse-amplitude modulation**

Class to generate and detect PAM signals.

Constructor of `M`-PAM signal with bit rate `Rb` and level spacing `level_spacing`. `pulse_shape` is generate using function [select_pulse_shape.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/signal/select_pulse_shape.m)
```
Mpam = PAM(M, Rb, level_spacing, pulse_shape);
```

Generating and detectin PAM symbols
```
symbols = Mpam.mod(dataTX);
dataRX = Mpam.demod(symbols);
```

Generating PAM signals with specified pulse shape
```
[signal, symbols] = Mpam.signal(dataTX);
```

- [QAM.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/signal/QAM.m): **Quadrature-Amplitude Modulation**

Class to generate and detect QAM signals. 

Constructor of `M`-QAM signal with bit rate `Rb`. `pulse_shape` is generate using function [select_pulse_shape.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/signal/select_pulse_shape.m)
```
Qam = QAM(M, Rb, pulse_shape);
```

Generating and detectin QAM symbols
```
symbols = Qam.mod(dataTX);
dataRX = Qam.demod(symbols);
```

`dataTX` can have multiple rows resulting in multidimensional symbols. For instance, for a dual-polarization signal `dataTX` would be a `2 x N` matrix.

Generating QAM signals with specified pulse shape
```
[signal, symbols] = Qam.signal(dataTX);
```

- [DPSK.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/signal/DPSK.m): **Differential Pulse-Shifting Keying**

Class to generate and detect DPSK signals. DPSK is a daugher class of class [QAM](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/signal/QAM.m)

Constructor of `M`-DPSK signal with bit rate `Rb`. `pulse_shape` is generate using function [select_pulse_shape.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/signal/select_pulse_shape.m)
```
Dpsk = DPSK(M, Rb, pulse_shape);
```

Generating and detectin DPSK symbols
```
symbols = Dpsk.mod(dataTX);
dataRX = Dpsk.demod(symbols);
```

`dataTX` can have multiple rows resulting in multidimensional symbols. For instance, for a dual-polarization signal `dataTX` would be a `2 x N` matrix.

Generating DPSK signals with specified pulse shape
```
[signal, symbols] = Dpsk.signal(dataTX);
```

- [OFDM.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/signal/OFDM.m): **Orthogonal frequency division multiplexing**

OFDM class can be used to generate double-side band (DSB), single-side band (SSB), DC-biased, and asymmetrically clipped (ACO) OFDM signals.

Contructor for OFDM with `Nc` subcarriers, `Nu` used subcarriers, nominal constellation size `CS`, bit rate `Rb`, prefix (either DSB, SSB, DC, ACO), and `powerAllocationType` (default 'Levin-Campello-MA')
```
Ofdm = OFDM(Nc, Nu, CS, Rb, prefix, powerAllocationType);
```

Example of generation and detection of OFDM signals using equalizer with parameters given in `eq`
```
[x, symbsTXm] = Ofdm.signal(Nsymb);
[Xn, AGCn, W] = Ofdm.detect(x, eq)
```

- [select_pulse_shape.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/signal/select_pulse_shape.m): **Select pulse shape**

Buils struct `pulse_shape` used to produce signals using the classes above. Options are either 'rect' (rectangular), 'rrc' (root-raised cosine), or 'rc' (raised cosine).

Function call:
```
pulse_shape = select_pulse_shape('rect', sps);
pulse_shape = select_pulse_shape('rrc', sps, rolloff, span);
pulse_shape = select_pulse_shape('rc', sps, rolloff, span);
```
`sps` denotes samples per symbol, 0 < rolloff < 1 is the roll-off factor, and `span` is the number of symbols over which the pulse shape spans.

- [Channels.m](https://github.com/Stanford-Optical-Comm-Group/optical-comm/blob/master/signal/Channels.m): **WDM channels**

Channels is a class of WDM channels. It contains information of wavelength, power, and propagation direction of all WDM channels. This class is used in submarine links optimization project. 

Constructor:
```
Ch = Channels(wavelength, power, direction)
```