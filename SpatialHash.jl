# SpatialHash.jl
# Provides a spatial hash for efficient neighbor lookup.

# A structure to hold the spatial hash
mutable struct SpatialHash
    cell_size::Float64
    grid::Dict{NTuple{3,Int}, Vector{Int}}  # Maps cell indices to lists of particle indices
end

# Constructor for SpatialHash
function SpatialHash(cell_size::Float64)
    SpatialHash(cell_size, Dict{NTuple{3,Int}, Vector{Int}}())
end

# Compute grid cell index for a given position
function cell_index(position::Vector{Float64}, cell_size::Float64)
    return (floor(Int, position[1] / cell_size),
            floor(Int, position[2] / cell_size),
            floor(Int, position[3] / cell_size))
end

# Insert a particle (by index) into the spatial hash
function insert!(sh::SpatialHash, particle_index::Int, particle::Particle)
    idx = cell_index(particle.position, sh.cell_size)
    if haskey(sh.grid, idx)
        push!(sh.grid[idx], particle_index)
    else
        sh.grid[idx] = [particle_index]
    end
end

# Retrieve neighbor particle indices by checking adjacent cells
function get_neighbors(sh::SpatialHash, particle::Particle)
    idx = cell_index(particle.position, sh.cell_size)
    neighbors = Int[]
    for dx in -1:1, dy in -1:1, dz in -1:1
        neighbor_cell = (idx[1] + dx, idx[2] + dy, idx[3] + dz)
        if haskey(sh.grid, neighbor_cell)
            append!(neighbors, sh.grid[neighbor_cell])
        end
    end
    return neighbors
end