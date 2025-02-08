# Particle-Based Viscoelastic Fluid Simulation in Julia

This repository contains a Julia implementation of a particle-based viscoelastic fluid simulation inspired by the work of Clavet, Beaudoin, and Poulin (2005):

> **Clavet, S., Beaudoin, P., & Poulin, P. (2005). Particle-based viscoelastic fluid simulation.**  
> In *Proceedings of the ACM SIGGRAPH/Eurographics symposium on Computer animation* (pp. 219–228).  
> DOI: [10.1145/1073368.1073400](https://doi.org/10.1145/1073368.1073400)

## Overview

This project implements a Lagrangian fluid simulation where the fluid is modeled as a set of interacting particles. Key techniques include:
- **Double Density Relaxation:** A method to enforce incompressibility by computing two types of local densities and applying corrective displacements.
- **Spring-Based Elasticity:** Particles are connected with springs to simulate viscoelastic behavior.
- **Spatial Hashing:** An efficient neighbor search is implemented using spatial hashing, which is essential for scalability.
- **Boundary Handling:** Fixed boundary collisions with sticky response to prevent particles from leaving the simulation domain.

## Project Structure

- **Constants.jl**  
  Contains simulation constants (e.g., interaction radii, stiffness parameters, visualization limits).

- **Particle.jl**  
  Defines the `Particle` type (with position, velocity, mass, and previous position) used in the simulation.

- **Spring.jl**  
  Defines the `Spring` type, which connects pairs of particles, along with parameters for stiffness, plasticity, and yield ratio.

- **SpatialHash.jl**  
  Implements a spatial hash data structure for efficient neighbor lookup.

- **FixedBoundaryCollisions.jl**  
  Contains functions to handle collisions with fixed boundaries, including a sticky collision response.

- **Update.jl**  
  Implements the main simulation update function following a prediction-relaxation scheme. This file applies gravity, viscosity, spring forces, double density relaxation, and boundary corrections. Finally, velocities are updated based on the displacement.

- **Visualize.jl**  
  Provides functions for visualizing the particles using Plots.jl.

- **Main.jl**  
  The main entry point for the simulation. It sets up the environment, instantiates particles (in this example, 100 particles with random positions), initializes springs and the spatial hash, and then runs the simulation for 500 frames to generate a GIF.



