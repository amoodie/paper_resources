This repository contains source code for *Optimized river diversion scenarios promote sustainability of urbanized deltas* (in press).  _Documentation for specific portions of the model code is included as line comments throughout the source code._

## Manuscript

Moodie, Andrew and Nittrouer, Jeffrey. *Optimized river diversion scenarios promote sustainability of urbanized deltas* (in press). Proceedings of the National Academy of Sciences.


## Delta model
The numerical delta model source code is located in [a sibling repository](https://github.com/amoodie/paper_resources/tree/master/Moodie_deltaiclobebuilding).
The model runs can be reproduced using the model source code and the script located in `model/scripts/`.


## Societal benefits and costs
The societal benefits and costs framework is implemented in the file `analysis/analysis_cost_benefit.m`
See main manuscript and supplementary material for derivation.

This script reproduces Figure 4 of the main text.

## Data

Model runs are too large to be preserved in any data repository (10--150 MB each, with 210 model runs).
This repository contains only summary data, which are located in `model/output_data`.

* `spinups/*` summary for spin-up simulations used as initial conditions for diversion simulations
* `engineered_single_analysis_summary.mat` contains a single table `eng` that summarizes all model simulation results.

**Note:** other files in this folder are redundant and included in the above summary file.

Within `engineered_single_analysis_summary.mat`, the following information is encoded:

* `eng.analysis` : cell array : containing summary information for each of 210 diversion simulations
* `eng.avulsion_location_list` : array : containing the normalized locations for each diversion simulation set (corresponds to rows of `eng.analysis`)
* `eng.Lblow` : scalar : backwater length estimate for YRd, Moodie et al., 2019 JGR
* `eng.tbl` : table : containing summary information for each set of diversion length simulation-replicate sets (i.e., mean and stddev of TA, LA, LL); data to reproduce Figure 2 in main text
* other arrays are non-aggregated datasets in this table; self-describing from variable name


## Acknowledgments and disclaimer
The model was created by Andrew J. Moodie as part of an National Science Foundation funded research project assessing the sustainability of anthropogenically influenced deltas.
The research is supported by Grant No. 1427262 and an NSF Graduate Research Fellowship to A.J.M. under Grant No. 1450681.
Any opinion, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
The author(s) guarantee no warranty or technical support for this model.
