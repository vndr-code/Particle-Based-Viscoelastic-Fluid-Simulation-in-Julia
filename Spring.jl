# Spring.jl
# Define a spring that connects two particles by their indices.

mutable struct Spring
    particle1::Int           # Index of the first particle
    particle2::Int           # Index of the second particle
    rest_length::Float64     # Current rest length of the spring
    stiffness::Float64       # Spring stiffness constant
    plasticity::Float64      # Plasticity constant (α)
    yield_ratio::Float64     # Yield ratio (γ)
end
