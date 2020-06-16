function [Z_rowcolval, V_rowcolval] = CalcZmatVvecLumpedFO(dof_data, num_dofs, mesh_data, interelem_VsrcZload)
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
% 2020-06-12: This function is edited for first order excitations. JT du
% Plessis

% Init:
num_lumped  = size(interelem_VsrcZload,1); % max (max possible number of lumped V's, max possible number of lumped Z's)
V_rowcolval = [ -ones(num_lumped,1)  -ones(num_lumped,1)  zeros(num_lumped,1) ];
Z_rowcolval = [ -ones(num_lumped,1)  -ones(num_lumped,1)  zeros(num_lumped,1) ];
Vnum        = 0; % actual number of V's (to be determined)
Znum        = 0; % actual number of Z's (to be determined)

% Cycle through all specified loads and sources and assemble their
% contributions:
for ii = 1:num_lumped

    % Find the global dof:
    candidate_dofs = dof_data.tri_to_dofs(interelem_VsrcZload(ii,1)).dofs; % select the first of the two tris arbitrarily, second would also work
    for jj = 1:size(candidate_dofs,2)
        this_dof = candidate_dofs(1,jj);
        if sort(dof_data.basis_supports(this_dof,:),2) == sort(interelem_VsrcZload(ii,1:2),2) % this is the dof
            
            % Edge length:
            this_edge   = dof_data.dofs_to_edges(this_dof,1);
            edge_verts  = mesh_data.edges(this_edge,1:2); % global nodes of current edge
            edgevectemp = mesh_data.node_coords(edge_verts(2),:) - mesh_data.node_coords(edge_verts(1),:);
            ELength     = sqrt(edgevectemp*edgevectemp');
            
            % Add V_vec contribution:
            V_val = interelem_VsrcZload(ii,3) + 1i*interelem_VsrcZload(ii,4);
            if abs(V_val) > 1e-8
                if dof_data.basis_supports(this_dof,1) == interelem_VsrcZload(ii,1)
                    dofsign =  1;
                else
                    dofsign = -1;
                end
                Vnum                = Vnum + 1;
                V_rowcolval(Vnum,:) = [ this_dof  1  dofsign*ELength*V_val ];
            end
            
            % Add Z_mat contribution:
            Z_val = interelem_VsrcZload(ii,5) + 1i*interelem_VsrcZload(ii,6);
            if abs(Z_val) > 1e-8
                Znum                = Znum + 1;
                Z_rowcolval(Znum,:) = [ this_dof  this_dof  (ELength^2)*Z_val ]; % lack of minus sign here is due to E_tot sign in EFIE: - E_scat + E_tot = E_inc 
%                 disp('');
            end
        end
    end




end