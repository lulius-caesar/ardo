Augmented Reality Design Optimization (ARDO)
============================================

This repository is for developing a reduced order model (ROM) for compsite
materials modeling with any given surface. It will try to optimize the laminate
configuration, so the number of plies, and the orienation and position of each
ply. The model will show in real time either the sugestion of draping process
or the consequences of draping at a particular orientation through an Augmented
Reality interphase. This allows the operator to determine the properties of the
properties of the part beforehand and make corrections in-situ, or ask for an
optimization of the laminate.

The program uses Classical Laminate Theory (CLT) to compute the strain and
stress on each node of a given 2D mesh with some albitrary boundary conditions. 

INPUT
-----

<!-- ARDO will use Leap motion for hand tracking -->
<!-- ADD ALSO the weibull modulus to compute the design allowables, considering
        a and b datasets. -->
- Mesh generated from Gmsh. It must include defined boundary conditions.
- Materials. These are given on the `materials.dat` file
- Laminate. These are given on the `laminate.dat` file.
  - [ ] The laminate can be optimized later with **ARDO** using *kinematic
        draping theory*
  - [ ] Generate laminate as chromosomes
- Boundary conditions. These are given on the `loads.dat` file. This file links
  the defined physical objects on the mesh file to apply displacements and loads
  (mechanical or hygro-thermal).


### Materials (materials.dat)
This file ignores new lines and `#`, so we can comment to make it more
understandable. The a material is defined in a single line with two functions.
1. Starting with `A` means that the material is *anisotropic* and needs to have
   the properties in the plane, i.e., 
   (E1 E2 G12 nu12 a1 a2 b1 b2 Xt Xc Yt Yc S12 S23).
2. Starting with `I` on the other hand, means that the material is isotropic and
   we just need to declare (E nu).


```
# Laminate configuration
1 0.127 0
0 0.254 90
0 0.254 90
1 0.127 0
```
```
10e2 0e2 0e0 0e0 10e0 0e0
```
```
# Material properties
# I = "isotropic"
# A = "anisotropic"

# Rubber
I 69e3 0.45

# Carbon
A 69e6 6e6 0.354 3e6 2.1e-3 2.1e-3 0.01 0.35 47e3 14e3 24e3 18e3 75e3 41e3 2e4

# Steel
I 210e3 0.354 
```
