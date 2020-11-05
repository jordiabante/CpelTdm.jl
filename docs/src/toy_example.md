# Toy Example

![CI](https://github.com/jordiabante/CpelTdm.jl/workflows/CI/badge.svg)
![Docs](https://github.com/jordiabante/CpelTdm.jl/workflows/Docs/badge.svg)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/jordiabante/CpelTdm.jl/blob/master/LICENSE.md)

The package includes a small toy example consisting of an unmatched group
comparison, which should take approximately 20 seconds to run. In particular, 
samples in group 1 were generated under an uncorrelated CPEL model with μ1=0.9,
while samples in group 2 were generated under an uncorrelated CPEL model with 
μ2=0.1. Thus, the example in particular consists of a region where the two groups 
strongly differ in mean methylation level (MML) and probability distribution of 
methylation (PDM), while they do not show any normalized methylation entropy 
(NME) difference. To run the toy example, go to the package folder and simply run 
`julia test/toy_example.jl`. CpelTdm will generate a total of 35 output files consisting 
of:

* A file per sample containing MML μ (16 total).
* A file per sample containing NME h (16 total).
* A file per differential analysis type (3 total).

in folder `test/cpeltdm/`.
