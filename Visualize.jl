# Visualize.jl
# Constants for visualization limits and marker size
const marker_scale = 6
const X_LIMITS = (-2.0, 2.0)
const Y_LIMITS = (-2.0, 2.0)
const Z_LIMITS = (0.0, 2.0)

# Function to visualize particles in 3D
function plot_particles(particles)
    p = scatter3d(
        [p.position[1] for p in particles],
        [p.position[2] for p in particles],
        [p.position[3] for p in particles],
        marker = :circle,
        markersize = marker_scale,
        xlims = X_LIMITS, ylims = Y_LIMITS, zlims = Z_LIMITS,
        aspect_ratio = :equal,
        legend = false,
        title = "Particles"
    )
    
end
