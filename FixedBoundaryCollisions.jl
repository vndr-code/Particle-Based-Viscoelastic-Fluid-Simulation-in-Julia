# FixedBoundaryCollisions.jl
#
# This file implements collision handling for a fixed (static) boundary.
# The boundary does not move (v_p = 0) so the collision response only affects the particle.
#
# For each particle that is outside the allowed limits, we:
# 1. Compute the current particle velocity (using current and previous positions).
# 2. Calculate the relative velocity (which is just the particle's velocity since the boundary is fixed).
# 3. Decompose that velocity into normal and tangential components.
# 4. Compute an impulse that cancels the normal component and reduces the tangential component using friction.
# 5. If the particle is very close to the boundary (penetration distance d < d_stick), compute an extra attractive "stickiness" impulse.
# 6. Update the particle's velocity by subtracting the total impulse.
# 7. Reset the particle's position to the boundary.
#
# Assumed constants (e.g., defined in Constants.jl):
#   const X_LIMITS = (0.0, 2.0)
#   const Y_LIMITS = (0.0, 2.0)
#   const Z_LIMITS = (0.0, 2.0)
#
# New constants for collision:
const d_stick = 0.1        # Threshold distance below which stickiness applies.
const k_stick = 50.0       # Strength of the stickiness (attraction).
const μ = 0.2              # Friction parameter (0 = full slip, 1 = no slip).

# Function to handle collisions of particles with fixed boundaries.
# Gravity is assumed to be along the z-axis (i.e. particles fall in z).
function resolveFixedBoundaryCollisions!(particles, dt)
    # Loop over all particles.
    for p in particles
        # First, compute the current velocity (v_i) from positions:
        # (Assuming p.position_prev exists and is updated each timestep.)
        v_i = (p.position .- p.position_prev) ./ dt

        # For each axis (1: x, 2: y, 3: z), check if the particle is outside the limits.
        for axis in 1:3
            # Choose the appropriate limits.
            lower_limit = axis == 1 ? X_LIMITS[1] : axis == 2 ? Y_LIMITS[1] : Z_LIMITS[1]
            upper_limit = axis == 1 ? X_LIMITS[2] : axis == 2 ? Y_LIMITS[2] : Z_LIMITS[2]
            
            # Check lower boundary.
            if p.position[axis] < lower_limit
                # Compute penetration distance d (positive value).
                d = lower_limit - p.position[axis]
                # Normal vector points into the domain.
                n = zeros(3)
                n[axis] = 1.0

                # Relative velocity v̄ = v_i - v_p, but here v_p = 0.
                v̄ = v_i

                # Compute normal component of v̄.
                v̄_normal = (dot(v̄, n)) * n
                # Tangential component is what remains.
                v̄_tangent = v̄ - v̄_normal

                # Compute collision impulse:
                # This impulse cancels the normal velocity and reduces tangential velocity by friction.
                I = v̄_normal - μ * v̄_tangent

                # Compute extra stickiness impulse if penetration is small.
                I_stick = zeros(3)
                if d < d_stick
                    I_stick = - dt * k_stick * d * (1 - d/d_stick) * n
                end

                # Total impulse to be applied.
                I_total = I + I_stick

                # Update the particle's velocity by subtracting the impulse.
                p.velocity -= I_total

                # Snap the particle's position to the boundary.
                p.position[axis] = lower_limit
            # Check upper boundary.
            elseif p.position[axis] > upper_limit
                d = p.position[axis] - upper_limit
                n = zeros(3)
                n[axis] = -1.0  # Normal points into the domain.
                
                v̄ = v_i
                v̄_normal = (dot(v̄, n)) * n
                v̄_tangent = v̄ - v̄_normal

                I = v̄_normal - μ * v̄_tangent

                I_stick = zeros(3)
                if d < d_stick
                    I_stick = - dt * k_stick * d * (1 - d/d_stick) * n
                end

                I_total = I + I_stick

                p.velocity -= I_total
                p.position[axis] = upper_limit
            end
        end
    end
end
