# Update.jl
#
# This file implements the main simulation step (Algorithm 1).
# Each step is laid out with comments that explain what it does.
#

# Apply gravity to each particle's velocity.
function applyGravity!(particles, dt, g)
    for p in particles
        p.velocity .+= dt .* g
    end
end

# Placeholder: Apply viscosity impulses between particle pairs (Section 5.3).
function applyViscosity!(particles, dt)
    N = length(particles)
    for i in 1:(N-1)
        for j in (i+1):N
            # Compute distance vector and its norm between particles i and j.
            r_vec = particles[j].position - particles[i].position
            r = norm(r_vec)
            if r < h   # Only consider neighbors within the interaction radius.
                q = r / h
                r_hat = r_vec / r  # Unit vector from i to j.
                # Project the velocity difference along the radial direction.
                u = dot(particles[i].velocity - particles[j].velocity, r_hat)
                if u > 0
                    # Compute the impulse magnitude.
                    I = dt * (1 - q) * (sigma * u + beta * u^2)
                    # Apply the impulses equally and oppositely.
                    particles[i].velocity .-= 0.5 * I * r_hat
                    particles[j].velocity .+= 0.5 * I * r_hat
                end
            end
        end
    end
end


# Save current positions as previous positions.
function savePreviousPositions!(particles)
    for p in particles
        p.position_prev .= copy(p.position)
    end
end

# Advance particles to their predicted positions using current velocities.
function advancePositions!(particles, dt)
    for p in particles
        p.position .+= dt .* p.velocity
    end
end

# Placeholder: Adjust springs (Section 5.2).
# adjustSprings!(particles, springs, dt)
function adjustSprings!(particles, springs, dt)
    N = length(particles)
    # Loop over each unique pair (i < j)
    for i in 1:(N-1)
        for j in (i+1):N
            p_i = particles[i]
            p_j = particles[j]
            r_vec = p_j.position - p_i.position
            r = norm(r_vec)
            q = r / h
            if q < 1.0
                # Check if a spring between i and j already exists.
                found = false
                for s in springs
                    if (s.particle1 == i && s.particle2 == j) || (s.particle1 == j && s.particle2 == i)
                        found = true
                        break
                    end
                end
                # If no spring exists, add one with rest length equal to h.
                if !found
                    push!(springs, Spring(i, j, h, spring_stiffness, plasticity, gamma_stretch))
                end
            end
        end
    end

    # Update rest lengths of existing springs based on deformation.
    for s in springs
        i, j = s.particle1, s.particle2
        p_i = particles[i]
        p_j = particles[j]
        r = norm(p_j.position - p_i.position)
        
        # Use gamma_stretch if the spring is stretched and gamma_compress if compressed.
        if r > s.rest_length
            d = gamma_stretch * s.rest_length  # Tolerable stretch deformation.
            if r > s.rest_length + d
                s.rest_length += dt * s.plasticity * (r - s.rest_length - d)
            end
        elseif r < s.rest_length
            d = gamma_compress * s.rest_length  # Tolerable compression deformation.
            if r < s.rest_length - d
                s.rest_length -= dt * s.plasticity * (s.rest_length - d - r)
            end
        end
    end

    # Remove springs whose rest length exceeds the interaction radius.
    delete_indices = []
    for (idx, s) in enumerate(springs)
        if s.rest_length > h
            push!(delete_indices, idx)
        end
    end
    for idx in reverse(delete_indices)
        deleteat!(springs, idx)
    end
end


# Placeholder: Apply spring displacements (Section 5.1).
function applySpringDisplacements!(particles, springs, dt)
    for s in springs
        # Get indices of the two connected particles.
        i = s.particle1
        j = s.particle2
        
        # Retrieve the particles.
        p_i = particles[i]
        p_j = particles[j]
        
        # Compute the vector and distance between them.
        r_vec = p_j.position - p_i.position
        r = norm(r_vec)
        if r == 0.0
            continue  # Avoid division by zero.
        end
        
        # Unit vector from particle i to j.
        r_hat = r_vec / r
        
        # Compute displacement:
        # D = dt^2 * k_spring * (1 - L/h) * (L - r) * r_hat
        D = (dt^2) * s.stiffness * (1 - s.rest_length / h) * (s.rest_length - r) * r_hat
        
        # Apply half of the displacement to each particle.
        p_i.position .-= 0.5 * D
        p_j.position .+= 0.5 * D
    end
end


# Placeholder: Double density relaxation (Section 4).
function doubleDensityRelaxation!(particles, shash, dt)
    for i in eachindex(particles)
        p_i = particles[i]
        density = 0.0
        near_density = 0.0
        neighbors = get_neighbors(shash, p_i)
        for j in neighbors
            if j == i
                continue
            end
            p_j = particles[j]
            r_vec = p_j.position - p_i.position
            r = norm(r_vec)
            if r < h
                q = r / h
                density += (1 - q)^2
                near_density += (1 - q)^3
            end
        end
        
        P = k * (density - rho0)
        P_near = k_near * near_density
        
        dx = zeros(3)
        for j in neighbors
            if j == i
                continue
            end
            p_j = particles[j]
            r_vec = p_j.position - p_i.position
            r = norm(r_vec)
            if r < h && r > 0
                q = r / h
                r_hat = r_vec / r
                D = (dt^2) * (P * (1 - q) + P_near * (1 - q)^2) * r_hat
                p_j.position .+= 0.5 * D
                dx .-= 0.5 * D
            end
        end
        p_i.position .+= dx
    end
end



# For now, we implement resolveCollisions as our fixed boundary collisions.
# (This corresponds to the collision part in Algorithm 1.)
function resolveCollisions!(particles, dt)
    resolveFixedBoundaryCollisions!(particles, dt)
end

# Main update function implementing Algorithm 1.
function update_simulation!(particles, springs, shash, dt, g)
    # 1. Apply gravity.                                             -DONE
    applyGravity!(particles, dt, g)
    
    # 2. Apply viscosity impulses.                                  -DONE
    applyViscosity!(particles, dt)
    
    # 3. Save current positions (for later velocity update).        -DONE
    savePreviousPositions!(particles)
    
    # 4. Advance to predicted positions.                            -DONE
    advancePositions!(particles, dt)
    
    # 5. Adjust springs (adding, removing, updating rest lengths).  -DONE
    adjustSprings!(particles, springs, dt)
    
    # 6. Apply spring displacements.                                -DONE
    applySpringDisplacements!(particles, springs, dt)
    
    # 7. Apply double density relaxation.
    doubleDensityRelaxation!(particles, shash, dt)
    
    # 8. Resolve collisions (using fixed boundary collisions here). -DONE (boundaries)
    resolveCollisions!(particles, dt)
    
    # 9. Update velocities based on new positions.
    for p in particles
        p.velocity .= (p.position .- p.position_prev) ./ dt
    end
end
