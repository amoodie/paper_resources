_NOTE: source code files from this project will be uploaded into this folder when the manuscript has been submitted._

This file provides documentation to the wrapper used to explore the model develop​ed in "Modeling deltaic lobe-building cycles and avulsions of the Yellow River delta, China" (Moodie et al., in prep).

<img src="./private/demo_lobe.png" alt="Demo image of lobe growth" width="600" align="middle">

# 1.0 Wrapper call syntax
A call to the wrapper looks like this:
```
[s] = virtualdelta_wrapper(preAvulThresh, preAvulTrigg, mouSwitch, QwSwitch, QwNum)
```

- `preAvulThresh` is the threshold amount of channel bed aggradation necessary for avulsion. It is normalized to the bankfull flow depth (i.e., fractional).
- `preAvulTrigg` is the amount of overbank flow necessary for avulsion. It is normalized to the bankfull fow depth (i.e., fractional).
- `mouSwitch` is the control for whether, and in what manner, lobe progradation occurs. See section "Input options" below for available options.
- `QwSwitch` is the control for the water discharge curve used in the model run. See section "Input options" below for available options.
- `QwNum` is the input "shape factor" that impacts the selected water discharge curve in the model. See section "Input options" below for available options.

# 2.0 Input options
## preAvulThresh
The amount of aggradation of the channel bed required for avulsion can be any numeric. It should be noted that the value is normalized to the bankfull flow depth, so a `preAvulThresh` of unity implies aggradation of one complete channel depth before avulsion. 


# 3.0 Disclaimer
The model was created by Andrew J. Moodie as part of an National Science Foundation funded research project assessing the sustainability of anthropogenically influenced deltas.
The research is supported by Grant No. 1427262 and an NSF Graduate Research Fellowship to A.J.M. under Grant No. 1450681.
Any opinion, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
The author(s) guarantee no warranty or technical support for this model.
