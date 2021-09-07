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
- [Carlomagno et al. 2020](https://www.spiedigitallibrary.org/journals/Journal-of-Astronomical-Telescopes-Instruments-and-Systems/volume-6/issue-3/035005/METIS-high-contrast-imaging-design-and-expected-performance/10.1117/1.JATIS.6.3.035005.full), METIS high-contrast imaging: design and expected performance
- [Bowens et al. 2021](https://www.aanda.org/articles/aa/full_html/2021/09/aa41109-21/aa41109-21.html), Exoplanets with ELT-METIS: I. Estimating the direct imaging exoplanet yield around stars within 6.5 parsecs

## Dependencies
HEEPS requires a [![Python 3.6](https://img.shields.io/badge/Python-3.6-brightgreen.svg)]() (or above) environment with the following libraries installed:   
```proper``` v3.2.5 --> [download](https://sourceforge.net/projects/proper-library/files/)  
```vip_hci``` v0.9.11  
```skimage```  v0.14.2  
```numpy``` v1.19.4  
```astropy```  v3.2.3  
```matplotlib```  v2.2.3  
```scipy```  v1.1.0  

You can use this tool to check all versions: [check_versions.ipynb](https://github.com/vortex-exoplanet/HEEPS/blob/master/tools/check_versions.ipynb)

## Quick start
This Jupyter Notebook will walk you through a simple HEEPS simulation: [demo.ipynb](https://github.com/vortex-exoplanet/HEEPS/blob/master/tools/demo.ipynb)