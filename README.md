# FluidFEM2
This is an implementation of the finite element method for the steady state, incompressible fluid flow.

### Method
This considers that the fluid's velocity can be described by a potential and that the curl of the velocity field is null. Leading to the equation:

$$
-\nabla^2 u = 0
$$

where $v = -\nabla u$

This requires robin boundary conditions, namely an "in" flux and the value of v at the "out" surface. Below is a simulation of a fluid moving around a "cylinder".

![fluidSimulation](https://github.com/user-attachments/assets/e75a66ef-492f-4b65-ae30-fccdb74b837b)
![fluidSimulation_HighRes](https://github.com/user-attachments/assets/6913bd90-fa99-452c-a9ec-2e176a1b85b7)

## Options
Users can set the geometry of the obstacle to either circles or rectangles, set their position and dimensions as well as the dimensions of the fluid container (the tube).
