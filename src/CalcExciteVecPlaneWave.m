function [V_vec] = CalcExciteVecPlaneWave(tri_dofs, node_coords, observer_map, k0, E_scalfac, theta_inc, phi_inc, eta_pol)
% The excitation vector is cacluated for the specfied dofs. <observer_map>
% is -1 for dofs not to be tested with, with ascending 1,2,3,... for the
% entries that has to be calculated. If only M entries must be calcuated,
% then <V_vec> only has M entries.
%
% The dof and mesh arguments are described here:
% mesh_data    : output of the function <CreateMeshData.m>
% dof_data     : output of the function <CreateBasisFunctions.m>
%
% 2019-12-16: Created. MMB.

% Init:
V_vec                = complex(zeros(size(observer_map,2),1));
% V_vec                = uint16(complex(size(ob
local_edge_nodes_def = [2 3
                        1 3
                        1 2];
numqpoints           = 7; % Assign the 7 point rule data: [Dunavant] degree of precision 5
qpoint1              = [ 0.333333333333333  0.333333333333333  0.333333333333333  0.225000000000000 ]; 
qpoint2              = [ 0.059715871789770  0.470142064105115  0.470142064105115  0.132394152788506 ]; 
qpoint3              = [ 0.797426985353087  0.101286507323456  0.101286507323456  0.125939180544827 ]; 
qruletri = zeros(7,4);
qruletri(1,1:4)      = [ qpoint1(1)  qpoint1(2)  qpoint1(3)  qpoint1(4) ];
qruletri(2,1:4)      = [ qpoint2(1)  qpoint2(2)  qpoint2(3)  qpoint2(4) ];
qruletri(3,1:4)      = [ qpoint2(2)  qpoint2(3)  qpoint2(1)  qpoint2(4) ]; 
qruletri(4,1:4)      = [ qpoint2(3)  qpoint2(1)  qpoint2(2)  qpoint2(4) ]; 
qruletri(5,1:4)      = [ qpoint3(1)  qpoint3(2)  qpoint3(3)  qpoint3(4) ]; 
qruletri(6,1:4)      = [ qpoint3(2)  qpoint3(3)  qpoint3(1)  qpoint3(4) ]; 
qruletri(7,1:4)      = [ qpoint3(3)  qpoint3(1)  qpoint3(2)  qpoint3(4) ];                     
unitvec_theta        =   [  cos(theta_inc)*cos(phi_inc)  cos(theta_inc)*sin(phi_inc)  -sin(theta_inc)  ];
unitvec_phi          =   [                -sin(phi_inc)                 cos(phi_inc)                0  ];
unitvec_beta         = - [  sin(theta_inc)*cos(phi_inc)  sin(theta_inc)*sin(phi_inc)   cos(theta_inc)  ];
unitvec_eta          = -cos(eta_pol)*unitvec_theta + sin(eta_pol)*unitvec_phi;


% Cycle over all tri_dofs (list of triangles, possibly listing some twice,
% with dofs associated with each listed triangle). Triangles are listed
% twice in case of junctions where two basis functions can be associated
% with the same edge of the same triangle. 
% Uver each triangle the relevant basis functions dot Einc are integrated
% and assembled into V-vec.
for ii = 1:size(tri_dofs,1)

    % Make a list of the indices into <V_vec> for the three basis functions
    % on this tri:
    triidx      = tri_dofs(ii,7:9);
    for jj = 1:3
        if triidx(1,jj) > 0 
            triidx(1,jj) = observer_map(1,triidx(1,jj));
        end
    end
    if max(triidx,2) < 0 % cycle the tri loop if no relevant basis to V_vec present in this triangle
        continue; 
    end
        
    % Calculate the contributions from this tri to <V_vec>:    
    nodtri       = tri_dofs(ii,1:3); % global nodes of current tri
    cross_val    = cross(node_coords(nodtri(2),:) - node_coords(nodtri(1),:), node_coords(nodtri(3),:) - node_coords(nodtri(1),:));
    TArea        = 0.5*sqrt(cross_val*cross_val');
    nodecoordmat = node_coords(nodtri,1:3);
    qcoords      = qruletri(:,1:3)*nodecoordmat;
    Einc_eval    = E_scalfac*exp(-1i*k0*qcoords*unitvec_beta')*unitvec_eta; % each row is the E-field(Ex,Ey,Ez) at each quad point
    for jj = 1:3 % cycle over the three edges of this triangle
        this_idx = triidx(1,jj);
        if this_idx > 0 % this basis function has a contribution to <V_vec>
            
            % Edge length:
            edge_verts  = tri_dofs(ii,local_edge_nodes_def(jj,1:2)); % global nodes of current edge
            edgevectemp = node_coords(edge_verts(2),:) - node_coords(edge_verts(1),:);
            ELength     = sqrt(edgevectemp*edgevectemp');
            
            % Eval basis function at the qpoints:
            RWG_eval = 0.5 * (ELength/TArea) * tri_dofs(ii,3+jj) * ( qcoords - ones(numqpoints,1)*node_coords(nodtri(1,jj),:) ); 
            
            % Calculate the integral result:
            int_eval = TArea * (  qruletri(:,4)'  ) * (  (Einc_eval.*RWG_eval)*[1 1 1]'  );
            
            % Assemble into <V_vec>:
            V_vec(this_idx) = V_vec(this_idx) + int_eval;
        end
    end
end
