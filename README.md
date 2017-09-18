# Two-Field Inflation Perturbation System

Mathematica package to numerically solve the two-field inflation system in both zero order and first order. It uses the parallel computing environment powered by MathLink to solve the systems of complex-valued partial differential equation that describe the scalar perturbations of two scalar fields and the metric tensor. It also solves the non-linear zero-order background system, which is fed into the calculation of the perturbations.

## Getting Started

Main560.m is the driver script. Change the directory paths of the source code and output data in this file (one needs to create the directory to store the output data in advance). The parameter scan range is also set here.

To run this script in Mathematica, in the front end, type

```
SetDirectory["/path/to/source/code"]
<< Main560.m
```

It will load in the packages: solveHybridPerturbation203.m, findAHybrid21.m, and FUN30.m. It will then run two scripts, ZeroOrderPart560.m and FirstOrderPart560.m, to solve the background and perturbation systems.

The spectra of the curvature perturbations of the modes scanned will be exported to the outputFolder (set in Main560.m).

The number of kernels used for parallel computing is by default the available cores in the system. One can also specify the number of kernels in FirstOrderPart560.m in

```
kernelList = LaunchKernels[numberOfKernels];
```

### Installing

Simply put all .m files into a directory, and specify that directory in Main560.m. Also create a directory to store the output data, and specify it in Main560.m. Then follow the instructions in Getting Started and you are all set.

## Built With

* [Mathematica](https://www.wolfram.com/mathematica/) - The Wolfram Mathematica Language

## Authors

* **Yu-Hsiang Lin** - [LinkedIn](https://www.linkedin.com/in/yhl2997925/)

## Publication

* Pisin Chen and Yu-Hsiang Lin, What initial condition of inflation would suppress the large-scale CMB spectrum? Phys. Rev. D 93, 023503 (2016). (First-author. Authors in alphabetical order.) [Phys. Rev. D](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.93.023503) [arXiv](https://arxiv.org/abs/1505.05980)
