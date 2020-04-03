function [quad_vertices_currents] = CalcElementsCurrentsQuad(I_vec, quad_dofs, quad_dofs_idx, num_quad, node_coords)
% This postprocessing routine calculates the values (Jx,Jy,Jz) at the three
% vertices of each triangle in the mesh, as evaluated from the perspective of that triangle.  
%
% Structure of <triangles_vertices_currents>:
% num_tri x 3(i.e. num vert per tri) x 3(three Cartesian components of J)
%
% 2019-12-16: Created. MMB.

% Init:
quad_vertices_currents = zeros(num_quad,4,3);
% local_edge_nodes_def        = [2 3
%                                1 3
%                                1 2];
% local_edge_nodes_def = [1 2; 2 3; 4 3; 1 4];
local_edge_nodes_def = [1 2;1 4 ; 4 3; 2 3];

% Cycle over all tri_dofs (list of triangles, possibly listing some twice,
% with dofs associated with each listed triangle). Triangles are listed
% twice in case of junctions where two basis functions can be associated
% with the same edge of the same triangle. 
for ii = 1:size(quad_dofs,1)
    this_quad  = quad_dofs_idx(ii,1); % triangle where this cycle's current contributions must be added
    n1        = quad_dofs(ii,1); % global nodes of current tri
    n2        = quad_dofs(ii,2); % global nodes of current tri
    n3        = quad_dofs(ii,3); % global nodes of current tri
    n4        = quad_dofs(ii,4);
    cross_val = cross(node_coords(n2,:) - node_coords(n1,:), node_coords(n3,:) - node_coords(n1,:));
    TArea     = 0.5*sqrt(cross_val*cross_val');
    
    b10 = -node_coords(n1,:) +node_coords(n2,:) - node_coords(n4,:) + node_coords(n3,:);
    b01 = -node_coords(n1,:) -node_coords(n2,:) + node_coords(n4,:) + node_coords(n3,:);
    drdv = 0.5 * b01;
    drdu = 0.5 * b10;
    J = cross(drdv,drdu);
    J_det = norm(J);
    n = J/J_det;
%     n_cross_l = cross(n,drdv);
%     n_cross_l1 = cross(n,drdu);
    
    for jj = 1:4 % cycle over the four edges of this quadrilateral
        this_dof = quad_dofs(ii,13-jj);
                
        if this_dof > 0 % a dof is associated with this edge and its contribution must be added to the vertex currents of this triangle
            edge_verts  = quad_dofs(ii,local_edge_nodes_def(jj,1:2)); % global nodes of current edge
%             if any(quad_dofs(ii,jj) ~= this_dof) && (quad_dofs(ii,jj) ~= -1) % If vertex is not part of the dof edge
%                t = 1;  
%             end
%             opp_indx = jj+2;
%             if opp_indx > 4
%                 opp_indx = 1;
%             end
%             opposite_edge = quad_dofs(ii,local_edge_nodes_def(opp_indx,:));
%             edge2 = quad_dofs(ii,mod(local_edge_nodes_def(jj,2)+1,4));
%             if edge1 == edge_verts(2)
%                 edge1 = quad_dofs(ii,mod(local_edge_nodes_def(jj,1)+1,4));
%             end
%             elen(1,:,1) = node_coords(edge_verts(1),:) - node_coords(opposite_edge(1),:);
%             elen(1,:,2) = node_coords(edge_verts(2),:) - node_coords(opposite_edge(2),:);
%             
%             
%             edgevectemp = node_coords(edge_verts(2),:) - node_coords(edge_verts(1),:);
%             ELength = norm(edgevectemp);
%             edgevectemp1 = edgevectemp/ELength;
%             ELength1 = norm(edgevectemp1);
%             Zeta = node_coords(edge_verts(1),:) - node_coords(quad_dofs(ii,jj),:);
%             ELength1 = 
%             ELength     = sqrt(edgevectemp*edgevectemp'); % TODO: get this to normalized coords
            
            switch jj
                case 1
%                     Zeta = norm(elen2) ;
                    zu = 1;zv = 0;
                    n_cross_l = cross(n,drdv);
                case 2
%                     Zeta = norm(elen1);
                    zv = 1;zu = 0;
                    n_cross_l = cross(n,drdu);
                case 3
%                     Zeta = norm(elen2);
                    zu = -1;zv = 0;
                    n_cross_l = cross(n,drdv);
                case 4
%                     Zeta = norm(elen1);
                    zu = 0;zv = -1;
                    n_cross_l = cross(n,drdu);
            end
            N = norm(n_cross_l);
            % Add to the vertex current totals for this triangle:
            for kk = 1:2 % contributions at the edge vertices, not the origin vertex 
                         % (note that basis direction must be incorporated via 
                         % appropriate tri_dofs data) 

                rho(1,1,1:3) =  ((zu*drdv) + (zv *drdu));
%                 Zeta = norm(elen(:,:,kk));
%                 rho(1,1,1:3) = node_coords(edge_verts(kk),:) - node_coords(quad_dofs(ii,jj),:);
%                 rho_edgevert(1,1,1:3) = node_coords(edge_verts(kk),:) - node_coords(quad_dofs(ii,jj),:); % origins of the three RWGs are the three tri nodes, in order
                quad_vertices_currents(this_quad,local_edge_nodes_def(jj,kk),1:3) = ... 
                    quad_vertices_currents(this_quad,local_edge_nodes_def(jj,kk),1:3) + ... 
                    I_vec(this_dof)*1* (N/J_det) * quad_dofs(ii,9-jj) * rho;
%                     I_vec(this_dof) * 0.5 * (ELength/TArea) * quad_dofs(ii,4+jj) * rho_edgevert;
            end
        end
    end
end
