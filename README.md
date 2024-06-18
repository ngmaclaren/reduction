This repository contains code to select optimized sentinel node sets which can closely approximate the state of all network nodes, as well as simulation and analysis code which produced the findings in our manuscript.

## optNS 

The `optNS` package is a small R package which implements our optimization algorithm; the name of the package is simply short for "optimize node sets." The `./optNS/` directory contains the code, which you can build with
```sh
R CMD build optNS
```
The built package is also included as `optNS_0.1-1.tar.gz`. To install, use
```sh
R CMD INSTALL optNS_0.1-1.tar.gz
```

The documentation for this package is a work in progress. Please contact me with any questions.

## Shell scripts for high-performance computing clusters

Our method should run in reasonable time for small networks on a standard laptop (I use a laptop with four Intel i3-5010U CPUs at 2.10GHz). In order to support larger networks we used a high-performance computing cluster to simulate node states across ranges of control parameters (our "ground truth" simulations) and to optimize node sets. The shell scripts we used to do this are in `./shell/` and will need to be modified for your parallel computing environment. 

## Dependencies

In addition to the `optNS` package, we use an in-house differential equation simulation package called [sdn](https://github.com/ngmaclaren/sdn). We use the following standard R repository packages:
- [deSolve](https://cran.r-project.org/package=deSolve)
- [ROI](https://cran.r-project.org/package=ROI)
- [ROI.plugin.qpoases](https://cran.r-project.org/package=ROI.plugin.qpoases)
- [igraph](https://cran.r-project.org/package=igraph)
- [parallel](https://cran.r-project.org/doc/manuals/r-release/fullrefman.pdf)
- [optparse](https://cran.r-project.org/package=optparse)
- [latex2exp](https://cran.r-project.org/package=latex2exp) (for the plots in `./analysis`)
