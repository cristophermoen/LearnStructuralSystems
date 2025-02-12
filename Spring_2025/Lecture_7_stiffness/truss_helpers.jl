function calculate_element_local_stiffness_matrix(E, A, L)

    ke_local =  zeros(Float64, (4, 4))
    ke_local[1, 1] = E * A / L 
    ke_local[1, 3] = -E * A / L
    ke_local[3, 1] = -E * A / L
    ke_local[3, 3] = E * A / L

    return ke_local

end


function define_rotation_matrix(θ)

    β = zeros(Float64, 4, 4)

    unit_rotation = [cos(θ) sin(θ)
                     -sin(θ) cos(θ)]

    β[1:2, 1:2] = unit_rotation
    β[3:4, 3:4] = unit_rotation

    return β

end


function calculate_element_length(node_i_coordinates, node_j_coordinates)

    L = norm(node_j_coordinates - node_i_coordinates)

    return L

end

function calculate_element_orientation(node_i_coordinates, node_j_coordinates)

    Δ = node_j_coordinates - node_i_coordinates
    
    θ = atan(Δ[2], Δ[1])

    return θ

end