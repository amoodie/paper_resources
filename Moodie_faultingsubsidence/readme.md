
## Manuscript

Moodie, Andrew and Paola Passalacqua. *When does faulting-induced subsidence drive distributary network reorganization?* Geophysical Research Letters. (2021). 


## Delta model
The core numerical delta model is implemented as *pyDeltaRCM* and housed in a [separate GitHub repository](https://github.com/DeltaRCM/pyDeltaRCM).
After installing *pyDeltaRCM*, the model runs can be executed with the scripts in this repository, located under `model/scripts`.

NOTE: to reproduce exact runs, use seed configurations in log files in `data/`.


## Data

* `data/Set_1/` all log files from simulation used in Set 1 simulations (Selenga-like configuration)
* `data/Set_2/` all log files from simulation used in Set 1 simulations (Mississippi-like configuration)
* `data/shapefiles` shapefiles used for volumetric estimates of sediment deposit in Selenga River delta, and making plots in Figure 1 of main text
* `data/realworld_measurements.ods` spreadsheet with measurements and compiled data to estimate displacement magnitude for real-world delta systems


## Acknowledgments and disclaimer

The research was supported by the National Science Foundation, by way of a Postdoctoral Fellowship to A.M. (EAR 1952772) and a grant to P.P. (EAR 1719670).
Any opinion, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
The author(s) guarantee no warranty or technical support for the model or analysis codes.
