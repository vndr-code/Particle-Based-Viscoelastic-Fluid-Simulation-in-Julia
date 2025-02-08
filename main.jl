#= 
Based on:
Clavet, Simon & Beaudoin, Philippe & Poulin, Pierre. (2005). 
Particle-based viscoelastic fluid simulation. 219-228. 
10.1145/1073368.1073400. 

Particle-Based Fluid Simulation:
It’s a method where the fluid is represented by particles 
(a Lagrangian approach), rather than a grid (Eulerian).

Double Density Relaxation:
Instead of using standard SPH kernels to compute forces, 
this method computes two kinds of “density” at each particle—one 
that considers all neighbors and another that focuses on very close neighbors. 
These densities are then used to adjust the particles’ 
positions so that the fluid remains incompressible and has a 
smooth surface with natural surface tension effects.

Comparison to SPH:
While it’s related to SPH in that both are particle-based, 
double density relaxation is a simpler and more robust alternative. 
It avoids complex curvature computations and instability issues often 
seen in traditional SPH.



 =#

 
using Pkg

# Activate environment
Pkg.activate(".")
Pkg.instantiate()
Pkg.add("Plots")
Pkg.add("StaticArrays")
Pkg.add("LinearAlgebra")

# main.jl
include("Constants.jl")
include("Particle.jl")
include("Spring.jl")
include("SpatialHash.jl")
include("Visualize.jl")
include("FixedBoundaryCollisions.jl")
include("SpatialHash.jl")
include("Update.jl")

using Plots
using LinearAlgebra
using Random

# Simulation parameters
dt = 0.01
g = [0.0, 0.0, -9.81]  # Gravity along the z-axis

# Create an array to store particles.
particles = Particle[]

# Create 100 particles with random positions within the limits:
# X in (-2,2), Y in (-2,2), Z in (1,2) so that they are above the floor.
for i in 1:100
    x = rand() * (2.0 - (-2.0)) + (-2.0)
    y = rand() * (2.0 - (-2.0)) + (-2.0)
    z = rand() * (2.0 - 1.0) + 1.0
    push!(particles, Particle([x, y, z], [0.0, 0.0, 0.0], 1.0))
end

# Create an empty springs array.
springs = Spring[]

# Create a spatial hash with cell size 1.0.
shash = SpatialHash(1.0)

# Define a helper function to update the spatial hash.
function updateSpatialHash!(sh::SpatialHash, particles)
    empty!(sh.grid)  # Clear the current grid.
    for (i, p) in enumerate(particles)
        insert!(sh, i, p)
    end
end

# Animation loop: each frame, update the simulation, update the spatial hash, then plot.
anim = @animate for frame in 1:500
    # Update simulation according to Algorithm 1.
    update_simulation!(particles, springs, shash, dt, g)
    
    # Update the spatial hash with new particle positions.
    updateSpatialHash!(shash, particles)
    
    # Plot the current particle positions.
    p = plot_particles(particles)
    title!(p, "Frame: $frame")
    p  # Return the plot for this frame.
end

gif(anim, "particle_simulation.gif", fps=30)

