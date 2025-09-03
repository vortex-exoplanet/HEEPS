# HEEPS
## The HCI End-to-End Performance Simulator
```python
				    )               (      (     
				 ( /(               )\ )   )\ )  
				 )\())  (     (    (()/(  (()/(  
				((_)\   )\    )\    /(_))  /(_)) 
				 _((_) ((_)  ((_)  (_))   (_))   
				| || | | __| | __| | _ \  / __|  
				| __ | | _|  | _|  |  _/  \__ \  
				|_||_| |___| |___| |_|    |___/
```
HEEPS is a high-contrast imaging (HCI) simulator, mostly geared towards studying the HCI performance of the ELT/[METIS](https://elt.eso.org/instrument/METIS/) thermal infrared instrument.

Currently, the simulator includes four coronagraphs:
- Classical Vortex Coronagraph (CVC)
- Ring Apodized Vortex Coronagraph (RAVC)
- Apodizing phase plate (APP)
- Classical Lyot Coronagraph (CLC)

## References
- [Delacroix, Absil, Orban de Xivry, et al. 2022](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/12187/121870F/The-High-contrast-End-to-End-Performance-Simulator-HEEPS/10.1117/12.2630341.short), The High-contrast End-to-End Performance Simulator (HEEPS): influence of ELT/METIS instrumental effects
- [Absil, Delacroix, Orban de Xivry, et al. 2022](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/12185/1218511/Impact-of-water-vapor-seeing-on-mid-infrared-high-contrast/10.1117/12.2627972.short), Impact of water vapor seeing on mid-infrared high-contrast imaging at ELT scale
- [Shinde, Delacroix, Orban de Xivry, et al. 2022](https://nanolithography.spiedigitallibrary.org/conference-proceedings-of-spie/12187/121870E/Modeling-the-vortex-center-glow-in-the-ELT-METIS-vortex/10.1117/12.2629855.short), Modeling the vortex center glow in the ELT/METIS vortex coronagraph
- [Bowens, Meyer, Delacroix, et al. 2021](https://www.aanda.org/articles/aa/full_html/2021/09/aa41109-21/aa41109-21.html), Exoplanets with ELT-METIS: I. Estimating the direct imaging exoplanet yield around stars within 6.5 parsecs
- [Carlomagno, Delacroix, Absil, et al. 2020](https://www.spiedigitallibrary.org/journals/Journal-of-Astronomical-Telescopes-Instruments-and-Systems/volume-6/issue-3/035005/METIS-high-contrast-imaging-design-and-expected-performance/10.1117/1.JATIS.6.3.035005.full), METIS high-contrast imaging: design and expected performance

## Dependencies
HEEPS depends on existing packages from the Python ecosystem, such as ``numpy``, ``scipy``, ``matplotlib``, ``astropy``, ``photutils``, ``skimage``.
It also relies on more specific package: [VIP](https://github.com/vortex-exoplanet/VIP) for HCI image processing and [PROPER](https://sourceforge.net/projects/proper-library/files/) for optical propagation

The versions under which the code is currently maintained are

[![Python](https://img.shields.io/badge/Python-3.12.7-brightgreen.svg)]()  
[![Proper](https://img.shields.io/badge/Proper-3.3.3-brightgreen.svg)]()
[![Scopesim](https://img.shields.io/badge/Scopesim-0.10.0-brightgreen.svg)]()  
[![Numpy](https://img.shields.io/badge/Numpy-1.26.4-brightgreen.svg)]()
[![Scipy](https://img.shields.io/badge/Scipy-1.15.0-brightgreen.svg)]()
[![Matplotlib](https://img.shields.io/badge/Matplotlib-3.10.6-brightgreen.svg)]()  
[![Vip-hci](https://img.shields.io/badge/Vip_hci-1.6.2-brightgreen.svg)]()
[![Photutils](https://img.shields.io/badge/Photutils-2.1.0-brightgreen.svg)]()
[![Astropy](https://img.shields.io/badge/Astropy-6.1.7-brightgreen.svg)]()
[![Skimage](https://img.shields.io/badge/Skimage-0.24.0-brightgreen.svg)]()

You can use a HEEPS tool to check all corresponding versions: [check_versions.ipynb](https://github.com/vortex-exoplanet/HEEPS/blob/master/notebooks/check_versions.ipynb)

## Quick start
This Jupyter Notebook will walk you through a simple HEEPS simulation: [demo.ipynb](https://github.com/vortex-exoplanet/HEEPS/blob/master/notebooks/demo.ipynb)
