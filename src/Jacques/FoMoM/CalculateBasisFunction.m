function [Rho] = CalculateBasisFunction(mesh_data, reduced_tri_dofs, edge, tri_num)

% 'Rho' is a 2x2 matrix containing both the RWG and Linear function magnitudes
% at both vertex points
% [RWG_v1 Lin_v1; RWG_v2 Lin_v2]

% Edge is local edge: 1,2,3
% Tri num is global number of triangle


Rho = zeros(2);

local_edge_nodes_def = [2 3
                        1 3
                        1 2];
% qpoint1              = [ 0.333333333333333  0.333333333333333  0.333333333333333  0.225000000000000 ]; 
% qpoint2              = [ 0.059715871789770  0.470142064105115  0.470142064105115  0.132394152788506 ]; 
% qpoint3              = [ 0.797426985353087  0.101286507323456  0.101286507323456  0.125939180544827 ]; 

for vertex = 1:2
    tri_point = mesh_data.node_coords(reduced_tri_dofs(tri_num,local_edge_nodes_def(edge,vertex)), :);
    mesh_TR = triangulation(mesh_data.tri_nodes,mesh_data.node_coords);
    Zeta = cartesianToBarycentric(mesh_TR,tri_num,tri_point); % Get local coord of vertex point
    
    
    for func_type = 1:2
        n1        = reduced_tri_dofs(tri_num,1); % global nodes of current tri
        n2        = reduced_tri_dofs(tri_num,2); % global nodes of current tri
        n3        = reduced_tri_dofs(tri_num,3); % global nodes of current tri
        cross_val = cross(mesh_data.node_coords(n2,:) - mesh_data.node_coords(n1,:), mesh_data.node_coords(n3,:) - mesh_data.node_coords(n1,:));
        
        t(1,:) = mesh_data.node_coords(n3,:) - mesh_data.node_coords(n2,:);
        t(2,:) = mesh_data.node_coords(n1,:) - mesh_data.node_coords(n3,:);
        t(3,:) = mesh_data.node_coords(n2,:) - mesh_data.node_coords(n1,:);
        
        TArea     = 0.5*sqrt(cross_val*cross_val');
        TLen      = norm(t(mod(edge+2,2)+1));
        
        t_index = mod(edge+1,3)+1;
        l_index = mod(edge,3)+1;
        first_order = -1;
        if func_type == 1
            sign = reduced_tri_dofs(tri_num,edge+3);
        else
            sign = 0.5*reduced_tri_dofs(tri_num,edge+9);
            first_order = 1;
        end
        
        % Only getting the magnitude of the basis function at the vertex
        Rho(vertex,func_type) = norm(sign*(TLen/(2*TArea))*(t(t_index,:)*Zeta(l_index)) + (first_order*(t(l_index,:)*Zeta(t_index))));
%         Rho(vertex,func_type) = sqrt(abs(sign*(TLen/(2*TArea))*(t(t_index,:)*Zeta(l_index)) + (first_order*(t(l_index,:)*Zeta(t_index)))));
        
    end

end
%     this_tri  = tri_dofs_idx(ii,1); % triangle where this cycle's current contributions must be added

            
end