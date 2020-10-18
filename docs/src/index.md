![CI](https://github.com/jordiabante/CpelTdm.jl/workflows/CI/badge.svg)
![Docs](https://github.com/jordiabante/CpelTdm.jl/workflows/Docs/badge.svg)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/jordiabante/CpelTdm.jl/blob/master/LICENSE.md)

## Description

CpelTdm is a julia package designed to perform DNA methylation differential analysis. 
The package is based on a recently published method [1], and has several advantages when 
compared to traditional DNA methylation differential analysis packages. CpelTdm not only 
detects mean methylation level differences between groups, but it can also detect 
significant differences in methylation entropy as well as in the probability
distribution of methylation (see [1] for technical details).

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

### Updating

`CpelTdm` and dependencies can be updated via the following command in julia's REPL:

```julia
(v1.3) pkg> update CpelTdm
```

## Running the tests

In order to test the package to ensure that it has been properly installed,
run the following command in a `julia` session:

```julia
(v1.3) pkg> test CpelTdm
```

If the package has been properly installed, then all tests will be successful.

## Authors

* **Jordi Abante**

## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/jordiabante/CpelTdm.jl/blob/master/LICENSE.md)
file for details.

## References

[1] Abante, J., Goutsias, J., CpelTdm.jl: a Julia package for targeted differential 
DNA methylation analysis, *bioArxiv* (2020), 343020.
