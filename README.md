# FluidFEM2

The mixed finite element method is implemented, solving the stokes equation for the benchmark lid-driven cavity problem. It is based on the Mats G. Larson. 2013 book "The finite element method: Theory, implementation and applications" (an adaptation of their version is also present in this repository for comparison).

The fluid is in steady state, in-compressible and viscous. The velocity field is zero at the walls and floor, and has "no-slip" boundary conditions with the lid which has a velocity [1,0].

Previously, this repository had a simulation of a fairly unrealistic fluid, considered as in-compressible and non-viscous, described by a potential "u". That code is still present in the directory "Irrotational_Fluids".


### Method
The mixed finite element method considers two different basis functions. Here, the pressure is approximated by piece-wise constant functions and the velocity field by the Crouzeix-Raviart function (which has an explicit form using the simple hat functions).

$$
-\nabla^2 \vec{u} + \nabla p = \rho \vec{f}
$$

where u is the velocity field and p is the pressure, f is a load force.

The lid-driven cavity problem has straight forward boundary conditions: the velocity field is zero on every surface except on the lid, where it is equal to the velocity of the lid.

### Lid-driven cavity problem
 
![velocity_pressure](https://github.com/user-attachments/assets/9488b18e-b8a2-4298-8e7d-f83da717a468)


### Non-viscous problem solution
![fluidSimulation](https://github.com/user-attachments/assets/e75a66ef-492f-4b65-ae30-fccdb74b837b)
![fluidSimulation_HighRes](https://github.com/user-attachments/assets/6913bd90-fa99-452c-a9ec-2e176a1b85b7)
