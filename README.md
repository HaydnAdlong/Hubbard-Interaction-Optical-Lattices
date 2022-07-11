# Microscopic calculation of Hubbard parameters for ultracold atoms for quantum gas microscopes
Calculation of the Hubbard on-site interaction U for cold atoms for quantum gas microscopes. This package is written in *Mathematica* and is currently applicable to quasi-1D and quasi-2D square optical lattices. For technical details on the code, please see the associated research paper freely available on the arXiv.
 

 ## Code units and notation
In this package we use units where `m=hbar=d=1` with
* `m` the particle mass,
* `hbar` the reduced Planck's constant and
* `d` the lattice spacing.

In these units the recoil energy is $`Vrec=1`$.

The standard notation includes
* `q` the centre of mass momentum
* `p` the incoming relative momentum
* `k` the relative momentum for integrating
* `tSigmaDirection` the hopping parameter of spin-sigma species along a direction (e.g. x or y in quasi-2D)



 ## Getting started
 To begin decide on the relevant package:
 - **1DHubbardParameters** for quasi-1D optical lattices
 - **2DHubbardParameters** for quasi-2D square optical lattices

Download the relevant package and then the package can then be loaded into a Mathematica notebook using, e.g.,
```
Get[<path-to-1DHubbardParameters.wl>];
```
 
 
 
 ## Citing this package
This package
