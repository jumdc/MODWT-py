# MODWT-py
Implementation of Maximum Overlap Discrete Wavelet Transform in Python

## Introduction 
This implementation follows the pyramid Algorithm as described in [1]

## Considerations
### Choice of filters
Generally the choice of filters need to be balance between two considerations : 
- wavelets with very short width can introduce undiserable artefacts into the resulting analyses. 
- wavelet filters with large L can be better however more coefficients are influcend by boundary effect, more computationally expensive.


## References
[1] Wavelet Methods for Time Series Analysis, Percival and Walden, 2000
