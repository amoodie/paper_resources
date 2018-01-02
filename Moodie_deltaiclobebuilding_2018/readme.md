NOTE: source code files from this project will be uploaded when the manuscript has been submitted.

This file provides documentation to the wrapper used to explore the model develop​ed in "Modeling deltaic lobe-building cycles and avulsions of the Yellow River delta, China" (Moodie et al., in prep).

[demo image](./private/demo_lobe.png)

# 1.0 Wrapper call syntax
A call to the wrapper looks like this:
```
[s] = virtualdelta_wrapper(preavulThresh, preavulTrigg, mouSwitch, Qwcall, Qwnum)
```


## 1.1 Windows (Windows 10)
* Visit the 'download' folder of this repository ([https://github.com/amoodie/research_outreach/tree/master/flooding_risk/download](./download)) and download the 64-bit Windows 10 FloodingModule installer.
* Run the installer by navigating to the download directory and executing the installer by double clicking. Note that you may have to run this as administrator, depending on where you intend to install the program and Matlab MCR.
* Follow the default options. If you must make any changes during the install process, be sure to note them (e.g., the install location), as you will may this information to run the module program.
* You can now run the program by navigating to the install location. 
  * note that the program may take up to a minute to launch the first time.

## 1.2 Mac (OSX)
forthcoming...

## 1.3 Linux (Ubuntu)
* Visit the 'download' folder of this repository ([https://github.com/amoodie/research_outreach/tree/master/flooding_risk/download](https://github.com/amoodie/research_outreach/tree/master/flooding_risk/download)) and download the 64-bit Linux FloodingModule installer.
* Run the installer by `cd`ing into the download location and running `./FloodingModule_lnx64.install`. Note that you may have to run this as root (i.e., with `sudo`) depending on where you intend to install the program and Matlab MCR.
* Follow the default options. If you must make any changes during the install process, be sure to note them (e.g., the install location), as you will need this information to run the module program.
* You can now run the program with `<FloodingModuleRoot>/application/run_FloodingModule.sh <MCRpath>`.
  * or, if you followed the default options, `/usr/local/FloodingModule/application/run_FloodingModule.sh /usr/local/MATLAB/MATLAB_Runtime/v92/`.
  * note that the program may take up to a minute to launch the first time.

## 1.4 Run from source
* The source Matlab code (`.m` file) is included in the repository in the 'source' folder. 
Two options for running the module
1. Save the file to a folder on your Matlab path and run `backwater_GUI()` from the Matlab console.
2. Open the `backwater_GUI.m` file in Matlab and select 'run'.


# 2.0 Module worksheet information

Worksheets are being written to accompany the GUI module. The aim of the worksheets is to help guide a discussion about flooding on deltas or rivers, as may be discussed in Earth Science courses. The modules currently available are targeted at specific age groups:

* High School Earth Science -- 9th to 12th graders

**Educators:** please let me know if you used this module, and if there are any ways you can see for us to help faciliate this sort of module in the future. Your feedback is very much appreciated. Email to Andrew Moodie, amoodie@rice.edu.


# 3.0 Disclaimer

This module utilizes a Matlab Compiled Runtime (MCR) program for backend calculations and plotting of the model results.
The program accompanying this module is to be run as a standalone application, but relies on the Matlab MCR; the program can be be downloaded at [http://www.coastalsustainability.rice.edu/outreach/](http://www.coastalsustainability.rice.edu/outreach/).
For help troubleshooting the Matlab executable, please visit [https://www.mathworks.com/products/compiler/mcr.html](https://www.mathworks.com/products/compiler/mcr.html).
The source code can be viewed at [https://github.com/amoodie/research_outreach](https://github.com/amoodie/research_outreach).

The module was created by Andrew J. Moodie as part of an National Science Foundation funded research project assessing the sustainability of anthropogenically influenced deltas.
The research is supported by Grant No. 1427262 and an NSF Graduate Research Fellowship to A.J.M. under Grant No. 1450681.
Any opinion, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
For more information, please visit [http://www.coastalsustainability.rice.edu/outreach/](http://www.coastalsustainability.rice.edu/outreach/).
