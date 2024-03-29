function [triangles_vertices_currents] = CalcElementsCurrentsFO(I_vec, tri_dofs, tri_dofs_idx, num_tri, node_coords, order, mesh_data)
% This postprocessing routine calculates the values (Jx,Jy,Jz) at the three
% vertices of each triangle in the mesh, as evaluated from the perspective of that triangle.
%
% Structure of <triangles_vertices_currents>:
% num_tri x 3(i.e. num vert per tri) x 3(three Cartesian components of J)
%
% 2019-12-16: Created. MMB.
% order = 0;
% Init:
triangles_vertices_currents = zeros(num_tri,3,3);
% local_edge_nodes_def        = [1 2
%                                1 3
%                                2 3];
local_edge_nodes_def        = [2 3
    1 3
    1 2];
if order == 1
    mesh_TR = triangulation(mesh_data.tri_nodes(1:2:end,:),mesh_data.node_coords);
else
    mesh_TR = triangulation(mesh_data.tri_nodes(:,:),mesh_data.node_coords);
end
% Somehow this function does not work as it should when not all the node
% coordinates are referenced by the triangulation

% w = warning('query','last');
% warning('off', w.identifier);
% if order == 1
%     tri_dofs_idx = tri_dofs_idx(1:2:end,1);
% end
% Cycle over all tri_dofs (list of triangles, possibly listing some twice,
% with dofs associated with each listed triangle). Triangles are listed
% twice in case of junctions where two basis functions can be associated

% with the same edge of the same triangle.

for ii = 1:size(tri_dofs,1)
    this_tri  = tri_dofs_idx(ii,1); % triangle where this cycle's current contributions must be added
    n1        = tri_dofs(ii,1); % global nodes of current tri
    n2        = tri_dofs(ii,2); % global nodes of current tri
    n3        = tri_dofs(ii,3); % global nodes of current tri
    cross_val = cross(node_coords(n2,:) - node_coords(n1,:), node_coords(n3,:) - node_coords(n1,:));
    TArea     = 0.5*sqrt(cross_val*cross_val');
    
    t(1,:) = node_coords(n3,:) - node_coords(n2,:);
    t(2,:) = node_coords(n1,:) - node_coords(n3,:);
    t(3,:) = node_coords(n2,:) - node_coords(n1,:);
    
    
    for jj = 1:(3+(3*order)) % cycle over the three edges of this triangle
        edge = mod(jj+1,3)+1;
%         edge = LocalEdgeNum(edge_nodes(ii,vertex),reduced_tri_dofs(ii,1:3));
        %         edge = jj;
        %         edge = jj;
        if jj > 3
            this_dof = tri_dofs(ii,12+edge);
            sign_iter = 6;
        else
            this_dof = tri_dofs(ii,6+edge);
            sign_iter = 0;
        end
        
        if this_dof > 0 % a dof is associated with this edge and its contribution must be added to the vertex currents of this triangle
            edge_verts  = tri_dofs(ii,local_edge_nodes_def(edge,1:2)); % global nodes of current edge
            edgevectemp = node_coords(edge_verts(2),:) - node_coords(edge_verts(1),:);
            ELength     = sqrt(edgevectemp*edgevectemp');
%                         ELength     = norm(t(mod(edge+2,2)+1));
%                         ELength     = norm(t(edge));
            first_order = -1;
            %             sign = tri_dofs(ii,3+sign_iter+(mod(jj,3)+1));
            sign = tri_dofs(ii,3+sign_iter+edge);
            if (abs(sign) == 2)
                first_order = 1;
                sign = 0.5*tri_dofs(ii,3+sign_iter+edge);
%                 ELength = 1;
%                 TArea = 1;
            end
            t_index = mod(edge+1,3)+1;
            l_index = mod(edge,3)+1;
            
            % Add to the vertex current totals for this triangle:
            for kk = 1:2 % contributions at the edge vertices, not the origin vertex
                % (note that basis direction must be incorporated via
                % appropriate tri_dofs data)
                 
                Zeta = cartesianToBarycentric(mesh_TR,this_tri,node_coords(edge_verts(kk),:)); % Get local coord of vertex point
                
                rho_edgevert(1,1,1:3) = (t(t_index,:)*Zeta(l_index)) + (first_order*(t(l_index,:)*Zeta(t_index)));
                triangles_vertices_currents(this_tri,local_edge_nodes_def(edge,kk),1:3) = ...
                    triangles_vertices_currents(this_tri,local_edge_nodes_def(edge,kk),1:3) + ...
                    I_vec(this_dof) * 0.5 * (ELength/TArea) * sign * rho_edgevert;
            end
        end
    end
end
