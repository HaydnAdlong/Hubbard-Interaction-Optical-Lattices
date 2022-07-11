# Microscopic calculation of Hubbard parameters for ultracold atoms for quantum gas microscopes
Calculation of the Hubbard on-site interaction U for cold atoms for quantum gas microscopes. This package is written in *Mathematica* and is currently applicable to quasi-1D and quasi-2D square optical lattices. For technical details on the code, please see the associated research paper freely available on the arXiv.
 
 ## Getting started
 To begin decide on the relevant package:
 - **1DHubbardParameters** for quasi-1D optical lattices
 - **2DHubbardParameters** for quasi-2D square optical lattices

Download the relevant package and then the package can then be loaded into a Mathematica notebook using, e.g.,
```
Get[<path-to-1DHubbardParameters.wl>];
```
 
 ## Code units and notation
In this package we use units where $`m=\hbar=d=1`$ with
* $`m`$ the particle mass,
* $`\hbar`$ the reduced Planck's constant and
* $`d`$ the lattice spacing.

In these units the recoil energy is $`V_{\text{r}}=1`$.

The standard notation includes
* $`q`$: centre of mass momentum
* $`p`$: incoming relative momentum
* $`k`$: relative momentum for integrating
* $`t_\sigma`$: hopping parameter of spin-\[Sigma] species
 
 
 
 ## Citing this package
This package
