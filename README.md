# Microscopic calculation of Hubbard parameters for ultracold atoms for quantum gas microscopes
Calculation of the Hubbard on-site interaction U for cold atoms for quantum gas microscopes. This package is written in *Mathematica* and is currently applicable to quasi-1D and quasi-2D square optical lattices. For technical details on the code, please see the associated research paper freely available on the arXiv.
 

 ## Code units and notation
In this package we use units where `m=hbar=d=1` with
* `m` the particle mass,
* `hbar` the reduced Planck's constant
* `d` the lattice spacing.

In these units the recoil energy is `Vrec=pi^2/2`.

The standard notation includes
* `q` the centre of mass momentum
* `p` the incoming relative momentum
* `k` the relative momentum for integrating
* `tSigmaDirection` the hopping parameter of spin-sigma species along a direction (e.g. x or y in quasi-2D)



 ## Using the package
 
 ### Getting started
 To begin decide on the relevant package:
 - **1DHubbardParameters** for quasi-1D optical lattices
 - **2DHubbardParameters** for quasi-2D square optical lattices

Download the relevant package and then the package can then be loaded into a Mathematica notebook using, e.g.,
```
Get[<path-to-1DHubbardParameters.wl>];
```

The packages enable the calculation of the Hubbard U term, as a function of 1D or 2D scattering length. This calculation depends on the optical lattice parameters, which in quasi-1D are:
* `vup` the depth of the optical lattice for the spin-up atoms
* `vdown` the depth of the optical lattice for the spin-down atoms
* `omegaperp` the trapping frequency of the 2D harmonic confinement.

The parameters in quasi-2D are:
* `vupx' the depth of the optical lattice for the spin-up atoms along the x dimension
* `vupy`, `vdownx` and `vdowny` follow the same definitions
* `omegaz` the trapping frequency of the 1D harmonic confinement.

As discussed in the paper the calculation is based on the exact calculation of the 2-particle scattering amplitude at (quasi) centre-of-mass momentum `q` and (quasi) relative on-shell momentum `pOnShell`. In quasi-2D, these parameters also involve a direction (i.e., `qx, qy, pOnShellx, pOnShelly`). Both parameters are input parameters in the code.

There are two options to performing this calculation depending on how much control the user wants. The key difference between these approaches is the control over the convergence parameters for the integrals and sums.

### Calculating U: The simple option
The simple option lumps all of the convergence parameters into a single parameter `convParameter`: the parameter takes on any non-zero integer, where a larger value yields a higher degree of convergence. Here is an example of using the code for quasi-1D
```
uFunc1D = setupHubbardUQuasi1D[vup = 12 Vrec, vdown = 12 Vrec, convParam = 4, pOnShell = 0.05, q = 0, omegaperp = 40]
```
Here, uFunc1D is a function of a 1D scattering length. For example, it can be plot using

### Calculating U: The involved option

 
 
 
 ## Citing this package
This package
