function calculate_element_local_stiffness_matrix(E, A, L)

    ke_local =  zeros(Float64, (4, 4))
    ke_local[1, 1] = E * A / L 
    ke_local[1, 3] = -E * A / L
    ke_local[3, 1] = -E * A / L
    ke_local[3, 3] = E * A / L

    return ke_local

end


function define_element_rotation_matrix(θ)

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


function define_global_element_stiffness_matrix(ke_local, β)

    ke_global = β' * ke_local * β

    return ke_global

end

function define_global_dof_for_element(node_number_i, node_number_j)

    dof_i = [node_number_i * 2 - 1, node_number_i * 2]
    dof_j = [node_number_j * 2 - 1, node_number_j * 2]

    dof = [dof_i; dof_j]   # [1, 2],  [3, 4]  => [1, 2, 3, 4]   

    return dof

end

function assemble_global_stiffness_matrix(num_dof, ke_global, element_global_dof)

    K = zeros(Float64, num_dof, num_dof)

    for i in eachindex(ke_global)

        K[element_global_dof[i], element_global_dof[i]] += ke_global[i]
        
    end

    return K

end
