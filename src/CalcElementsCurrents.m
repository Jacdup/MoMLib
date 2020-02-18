function [triangles_vertices_currents] = CalcElementsCurrents(I_vec, tri_dofs, tri_dofs_idx, num_tri, node_coords)
% This postprocessing routine calculates the values (Jx,Jy,Jz) at the three
% vertices of each triangle in the mesh, as evaluated fromthe perspective of that triangle.  
%
% Structure of <triangles_vertices_currents>:
% num_tri x 3(i.e. num vert per tri) x 3(three Cartesian components of J)
%
% 2019-12-16: Created. MMB.

% Init:
triangles_vertices_currents = zeros(num_tri,3,3);
local_edge_nodes_def        = [2 3
                               1 3
                               1 2];

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
    for jj = 1:3 % cycle over the three edges of this triangle
        this_dof = tri_dofs(ii,6+jj);
        if this_dof > 0 % a dof is associated with this edge and its contribution must be added to the vertex currents of this triangle
            edge_verts  = tri_dofs(ii,local_edge_nodes_def(jj,1:2)); % global nodes of current edge
            edgevectemp = node_coords(edge_verts(2),:) - node_coords(edge_verts(1),:);
            ELength     = sqrt(edgevectemp*edgevectemp');
            
            % Add to the vertex current totals for this triangle:
            for kk = 1:2 % contributions at the edge vertices, not the origin vertex 
                         % (note that basis direction must be incorporated via 
                         % appropriate tri_dofs data) 
                rho_edgevert(1,1,1:3) = node_coords(edge_verts(kk),:) - node_coords(tri_dofs(ii,jj),:); % origins of the three RWGs are the three tri nodes, in order
                triangles_vertices_currents(this_tri,local_edge_nodes_def(jj,kk),1:3) = ...
                    triangles_vertices_currents(this_tri,local_edge_nodes_def(jj,kk),1:3) + ... 
                    I_vec(this_dof) * 0.5 * (ELength/TArea) * tri_dofs(ii,3+jj) * rho_edgevert;
            end
        end
    end
end
