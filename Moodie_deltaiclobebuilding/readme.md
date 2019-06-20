**_NOTE: source code files from this project will be uploaded into this folder when the manuscript has been submitted._**

This file provides documentation to the wrapper used to explore the model develop​ed in "Modeling deltaic lobe-building cycles and avulsions of the Yellow River delta, China" (Moodie et al., in prep). See manuscript for a complete description of the model framework and parameterization. _Documentation for specific portions of the model code is included as line comments throughout the source code._

<img src="./private/demo_lobe.png" alt="Demo image of lobe growth" width="600" align="middle">

## Model wrapper
The model is split across many functions, organized together into categories below the `source` folder.
A call to the wrapper looks like this:

```
[s] = virtualdelta(cfg)
```

where `cfg` is the configuration structure, explained below.

A typical script will 1) set up the configuration structure `cfg`, and 2) call the model through the wrapper.
The bulk of the configuration structure is initialized with a call to a "system initializer", such as `[cfg] = yellowriver_configuration();`.
This function will return the structure with most of the required configurations specified, including delta slope, domain length, numerical stability parameters, sediment transport parameters, and depositional width parameters. 
To configure the model for a different system, you should create a copy of `yellowriver_configuration.m` and modify it under an appropriate name.

The script must then add a few "runtime" options to the configuration structure, which are the parameters that are varied in the model runs in the manuscript.
These values include, the superelevation threshold for avulsion and the water discharge scheme.
Of course, any configuration parameters initialized by the system initializer can be overwritten in the script at this point.


# 3.0 Acknowledgements and disclaimer
The model was created by Andrew J. Moodie as part of an National Science Foundation funded research project assessing the sustainability of anthropogenically influenced deltas.
The research is supported by Grant No. 1427262 and an NSF Graduate Research Fellowship to A.J.M. under Grant No. 1450681.
Any opinion, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
The author(s) guarantee no warranty or technical support for this model.
