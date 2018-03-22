#   Research on optical communication systems
> Developed during graduate school at Stanford University.

This project contains code for analyses and simulations of optical communications systems. 

## Folder description

- [components](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/components): components and devices.
- [docs](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/docs): documentation (github.io).
- [measurement](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/measurement): measurement and characterization.
- [processing](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/processing): signal processing and impairment mitigation.
- [projects](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects): projects developed using the base code.
- [propagate](hhttps://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/propagate): TBD.
- [signal](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/signal): classes on modulation formats.
- [tools](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/tools): TBD.
- [utilities](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/utilities): auxiliary functions.


# Projects description
- [examples](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/examples/): some simple examples on how to run simulations of basic optical communication systems

- [submarine_link_opt](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/submarine_link_opt/): optimization of submarine optical links limited by energy constraints. Simulations of erbium-doped fiber amplifiers. EDFAs are modeled using the Standard Confined-Doping (SCD) model and capacity optimization is performed using the particle swarm optimization algorithm.
  - [submarine_link_opt/doc/](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/submarine_link_opt/doc): documentation. Latex file containing some of the theoretical derivations and analyses
  - [submarine_link_opt/data/](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/submarine_link_opt/data): data for some erbium-doped fibers
  - [submarine_link_opt/f/](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/submarine_link_opt/f): auxiliary functions and classes used in submarine_link_opt/
  - [submarine_link_opt/results/](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/submarine_link_opt/results): folder reserved for saving files of simulations on cluster and posterior processing.
  - [submarine_link_opt/validation/](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/submarine_link_opt/validation): scripts used to test and validate functions in submarine_link_opt/

__Publications__

  *J. K. Perin, and J. Kahn, "Importance of Amplifier Physics in Maximizing the Capacity of Submarine Links," arXiv, 2018. [PDF](https://arxiv.org/abs/1803.07905)


- [dsp_free_coherent](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/dsp_free_coherent/): analysis and simulations of coherent and differentially coherent receivers. This includes DSP-based systems as well as systems based on analog signal processing
  - [dsp_free_coherent/analog](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/dsp_free_coherent/analog): functions and classes for modeling building blocks in the dsp_free_coherent receiver based on analog signal process
  - [dsp_free_coherent/analysis](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/dsp_free_coherent/analysis): analysis scripts. These are typically oversimplified simulations to better understand some concepts 
  - [dsp_free_coherent/doc](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/dsp_free_coherent/doc): documentation. Latex file containing some of the results and analysis
  - [dsp_free_coherent/f](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/dsp_free_coherent/f): auxiliary functions and classes used in dsp_free_coherent/
  - [dsp_free_coherent/results](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/dsp_free_coherent/results): folder reserved for saving files of simulations on cluster and posterior processing.
  - [dsp_free_coherent/validate](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/dsp_free_coherent/validate): scripts reserved for validation

__Publications__

  * J. Krause Perin, A. Shastri, J. Kahn, "DSP-Free Coherent Receivers for Data Center Links," OFC, 2018. [PDF](http://www.stanford.edu/~jkperin/OFC_DSP_free_coherent.pdf) 
  * J. Krause Perin, A. Shastri, J. Kahn, "Data Center Links Beyond 100 Git/s per Wavelength," Photonics West, 2018. [PDF](http://www.stanford.edu/~jkperin/PW_DC_review.pdf) 
  * J. Krause Perin, A. Shastri, J. Kahn, "Data Center Links Beyond 100 Git/s per Wavelength," Optical Fiber Technology, 2017. [PDF](http://www.stanford.edu/~jkperin/data_center_review.pdf) 
  * J. Krause Perin, A. Shastri, and J. Kahn, "Design of Low-Power DSP-Free dsp_free_coherent Receivers for Data Center Links," J. Lightw. Technol., vol. 35, no. 21, pp. 4650–4662, 2017. [PDF](http://www.stanford.edu/~jkperin/DSP-free_coherent.pdf)

- [apd/](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/apd): analysis and simulations of intensity-modulated direct-detected (IM-DD) optical systems using avalanche photodiodes  
  - [apd/doc](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/apd/doc): documentation. Latex file containing some of the theoretical derivations and analyses
  - [apd/f](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/apd/f): auxiliary functions used in apd/
  - [apd/results](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/apd/results): folder reserved for saving files of simulations on cluster and posterior processing.\

__Publications__

  * J. Krause Perin, M. Sharif, J.M. Kahn, "Sensitivity Improvement in 100 Gbit/s-per- Wavelength Links using Semiconductor Optical Amplifiers or Avalanche Photodiodes," J. Lightw. Technol., vol. 34, no. 33, pp. 5542–5553, 2016. [PDF](http://www.stanford.edu/~jkperin/SOA_vs_APD_100G.pdf)

- [ofdm/](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/ofdm): Comparison between single-carrier and multicarrier modulation formats for data center links
  - [ofdm/doc](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/ofdm/doc): documentation. Latex file containing some of the theoretical derivations and analyses
  - [ofdm/f/](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/ofdm/f): auxiliary functions used in ofdm/
  - [ofdm/analysis/](https://github.com/Stanford-Optical-Comm-Group/optical-comm/tree/master/projects/ofdm/analysis): analysis scripts. These are typically oversimplified simulations to better understand some concepts.

__Publications__

  * J. Krause Perin, M. Sharif, J. M. Kahn, "Modulation Schemes for Single-Wavelength 100 Gbits/s Links: Multicarrier," J. of Lightwave Technol., vol.33, no. 24, pp.5122-5132, Dec. 15, 2015. [PDF](http://ee.stanford.edu/~jmk/pubs/100.G.single-laser.multicarrier.JLT.15.pdf)
  * M. Sharif, J. Krause Perin, and J. M. Kahn, "Modulation Schemes for Single-Wavelength 100 Gbits/s Links: Single-Carrier," J. of Lightwave Technol., vol.33, no.20, pp.4268-4277, Oct. 15, 2015. [PDF](http://ee.stanford.edu/~jmk/pubs/100.G.single-laser.single-carrier.JLT.15.pdf)
