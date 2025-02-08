# Particle.jl

# Define a mutable Particle type with fields for position, velocity, mass, and position_prev.

mutable struct Particle
    position::Vector{Float64}      # Current position [x, y, z]
    velocity::Vector{Float64}      # Current velocity [vx, vy, vz]
    mass::Float64                  # Mass of the particle
    position_prev::Vector{Float64} # Previous position; initially same as position

    # Constructor: When creating a new Particle, position_prev is set to a copy of position.
    function Particle(position::Vector{Float64}, velocity::Vector{Float64}, mass::Float64)
        new(position, velocity, mass, copy(position))
    end
end

