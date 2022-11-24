# Microscopic calculation of Hubbard parameters for quantum gas microscopes
A package for the calculation of the Hubbard on-site interaction U for cold atoms in optical lattices. This package is written in *Mathematica* and is currently applicable to quasi-1D and (square) quasi-2D optical lattices. 

## How it works

The Hubbard U term is calculated such that the Hubbard scattering amplitude reproduces the exact scattering amplitude for two atoms in the limit of vanishing relative momentum. This packages automates the calculation and equating of the Hubbard and exact scattering amplitudes.

For more details please see the associated research paper freely available on the arXiv.

 ## Using the package
 
 ### Getting started
 To begin, decide on the relevant package:
 - **1DHubbardParameters** for quasi-1D optical lattices
 - **2DHubbardParameters** for quasi-2D (square) optical lattices

Download the relevant package and then the package can then be loaded into a Mathematica notebook using, e.g.,
```
Get[<path-to-1DHubbardParameters.wl>];
```

### Units
In this package the units are set by `m=hbar=d=1` with
* `m` the atom mass,
* `hbar` the reduced Planck's constant
* `d` the lattice spacing.


### Quasi-1D optical lattice
The key function in quasi-1D is `setupHubbardUQuasi1D[vup, vdown, convParam, pOnShell, omegaperp]` which outputs the Hubbard U as a function of `a1dinv` (the inverse 1D scattering length). The inputs of the function are:
* `vup`: depth of the optical lattice for the spin-up atoms [in units of `1/md^2`]
* `vdown`: depth of the optical lattice for the spin-down atoms [in units of `1/md^2`]
* `convParam`: convergence parameter for the integrals and sums, which takes on any non-zero integer (typically a value < 10 will suffice)
* `pOnShell`: relative quasi-momentum for equating the exact and Hubbard scattering amplitudes (must be taken as a limit to zero, but values around 0.05 typically suffice) [in units of `1/d`]
* `omegaperp` the trapping frequency of the 2D harmonic confinement [in units of `1/md^2`].

The following is an example of using the function, with experimentally feasible parameters:
```
uFunc1D = setupHubbardUQuasi1D[vup = 12 Vrec, vdown = 12 Vrec, pOnShell = 0.05, omegaperp = 40, convParam = 4]
```
The function `uFunc1D` can now be plotted
```
Plot[
 uFunc1D[a1dinv],
 {a1dinv, -30, 30},
 Frame -> True,
 FrameLabel -> {"d/a1D", "U/m d^2"}
 ]
```
yielding

<img width="360" alt="image" src="https://user-images.githubusercontent.com/93458010/203864975-0f0c101c-9f2b-4ceb-b8f8-89295492f93d.png">


One can also introduce the 3D and 1D scattering length relationship
```
lperp = 1/Sqrt[omegaperp];
a1dinvFunc[a3d_] := (-2 (lperp/2 (lperp/a3d + Zeta[1/2]/Sqrt[2])))^-1
```
Then, `uFunc1D` can now be plotted as a function of the 3D scattering length
```
Plot[
  uFunc1D[a1dinvFunc@a3d],
  {a3d, -1, 1},
  Frame -> True,
  FrameLabel -> {"a3d/d", "U/m d^2"}
  ]
```
<img width="360" alt="image" src="https://user-images.githubusercontent.com/93458010/203864943-3483c5e5-1a26-4e07-8805-d8bca6913af3.png">

The limit of zero relative quasi-momentum can be tested by reducing the value to `pOnShell` (e.g. `pOnShell = 0.025`), and the convergence can be tested by increasing the value of `convParam` (e.g. `convParam = 5`).




## old




In these units the recoil energy is `Vrec=pi^2/2`.

The standard notation includes
* `q` the centre of mass momentum
* `p` the incoming relative momentum
* `k` the relative momentum for integrating
* `tSigmaDirection` the hopping parameter of spin-sigma species along a direction (e.g. x or y in quasi-2D)




The packages enable the calculation of the Hubbard U term, as a function of 1D or 2D scattering length. This calculation depends on the optical lattice parameters, which in quasi-1D are:
* `vup` the depth of the optical lattice for the spin-up atoms
* `vdown` the depth of the optical lattice for the spin-down atoms
* `omegaperp` the trapping frequency of the 2D harmonic confinement.

The parameters in quasi-2D are:
* `vupx` the depth of the optical lattice for the spin-up atoms along the x dimension
* `vupy`, `vdownx` and `vdowny` follow the same definitions
* `omegaz` the trapping frequency of the 1D harmonic confinement.

As discussed in the paper the calculation is based on the exact calculation of the 2-particle scattering amplitude at (quasi) centre-of-mass momentum `q` and (quasi) relative on-shell momentum `pOnShell`. In quasi-2D, these parameters also involve a direction (i.e., `qx, qy, pOnShellx, pOnShelly`). The calculation requires the limit of `pOnShell` to zero.

There are two options to performing this calculation depending on how much control the user wants. The key difference between these approaches is the control over the convergence parameters for the integrals and sums.

 
 ## Citing this package
This package
