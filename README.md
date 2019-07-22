# HEEPS: High-contrast End-to-End Performance Simulator

					    )               (      (     
					 ( /(               )\ )   )\ )  
					 )\())  (     (    (()/(  (()/(  
					((_)\   )\    )\    /(_))  /(_)) 
					 _((_) ((_)  ((_)  (_))   (_))   
					| || | | __| | __| | _ \  / __|  
					| __ | | _|  | _|  |  _/  \__ \  
					|_||_| |___| |___| |_|    |___/    



## Overview
HEEPS is a high-contrast imaging (HCI) simulator, mostly geared towards studying the HCI performance of ELT instrument; "METIS". 

[![N|Solid](https://i2.wp.com/metis.strw.leidenuniv.nl/wp-content/uploads/2017/11/logo_with_text.png?resize=300%2C238)](https://i2.wp.com/metis.strw.leidenuniv.nl/wp-content/uploads/2017/11/logo_with_text.png?resize=300%2C238)

Currently simulator includes three coronagraphs:
- Classical vortex
- Ring apodized vortex coronagraph (RAVC)
- Apodizing phase plate (APP)

For technical document about the simulator see SPIE paper by Brunella et al: [link](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/9909/1/End-to-end-simulations-of-the-E-ELTMETIS-coronagraphs/10.1117/12.2233444.short)

## Getting started with HEEPS
1. Python 3.6 is required to run HEEPS package, if not already installed, we recommend downloading and installing Anaconda package for Python 3.6, see [link](https://www.anaconda.com/download/#linux). 
2. The 'PROPER' library can be downloaded from [link](https://drive.google.com/file/d/1CLsUjNI4Vfe8QdZ-qpEhc3DrgZb9H2_5/view?usp=sharing).
    (i) After downloading the "proper_v3.0d1_python_3.x_30jul18.zip" file, extract it.
    (ii) Go to the directory and execute "python setup.py install"  in a terminal
2. The 'VIP' library for ADI processing and generating contrast curves can be installed by executing "pip install vip_hci" in a terminal. Please see the GitHub [page]( https://github.com/vortex-exoplanet/VIP) for more details about the package.

## Using HEEPS

1. Download all the HEEPS file as a zip and extract it.
2. A default simulation parameters are described in a file called "read_config.py". These parameters can be modified/overwritten in a script.
3. See files "example_coronagraph_psf.py" and "ADI_processing.py" to get familiar with the module
4. See "manual.pdf" for more details about the architecture of the HEEPS package and its various sub-modules.
