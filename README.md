![CI](https://github.com/jordiabante/CpelTdm.jl/workflows/CI/badge.svg)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/jordiabante/CpelTdm.jl/blob/master/LICENSE.md)

## Description

CpelTdm is a julia package especifically desgined for haplotype allele-specific 
methylation based on the method in [1].

## Testing

CpelTdm is tested against Julia `1.3.0` on the latest versions of Linux, macOS and Windows.

## Getting Started

### Prerequisites

* julia v1.3.0
* git.

### Installing

`CpelTdm` and dependencies can be installed via the following command in julia's REPL:
```julia
(v1.3) pkg> add https://github.com/jordiabante/CpelTdm.jl.git
```

## Running the tests

In a `julia` session run
```julia
(v1.3) pkg> test CpelTdm
```

## Authors

* **Jordi Abante**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md)
file for details.

## References
[1] Abante, J., Goutsias, J., CpelTdm.jl: a Julia package for targeted diﬀerential DNA methylation analysis, *Bioinformatics* 2020 XYZ.
