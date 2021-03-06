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
Vnum        = -1; % actual number of V's (to be determined)
Znum        = -1; % actual number of Z's (to be determined)

% Cycle through all specified loads and sources and assemble their
% contributions:
for ii = 1:num_lumped
    
    if ii > 20
        test = 1;
    end

    % Find the global dof:
    candidate_dofs = dof_data.tri_to_dofs(interelem_VsrcZload(ii,1)).dofs; % select the first of the two tris arbitrarily, second would also work
    for jj = 1:size(candidate_dofs,2)
        this_dof = candidate_dofs(1,jj);
        
        
        % The following is a story of a little sorting algorithm
        % It goes something like this
        [~,b] = ismember(dof_data.basis_supports(this_dof,:),interelem_VsrcZload(:,1:2));
        if (nnz(b) == 2)
            [b(2),~] = ind2sub(size(interelem_VsrcZload(:,1:2)),b(2)); % Get actual index
            % Then, swap the elements of column 2 around to row of column 1
            temp = interelem_VsrcZload(b(1), 2);
            interelem_VsrcZload(b(1), 2) = interelem_VsrcZload(b(2),2);
            interelem_VsrcZload(b(2), 2) = temp;
            
        end
        
        % This check does not work if the rows are not aligned (if the
        % triangles do not match)
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
                Vnum                  = Vnum + 2;
                V_rowcolval(Vnum,:)   = [ this_dof-1   1  dofsign*ELength*V_val ];
                V_rowcolval(Vnum+1,:) = [ this_dof 1  0]; % Extra DOF for first order spec
            end
            
            % Add Z_mat contribution:
            Z_val = interelem_VsrcZload(ii,5) + 1i*interelem_VsrcZload(ii,6);
            if abs(Z_val) > 1e-8
                Znum                  = Znum + 2;
                Z_rowcolval(Znum,:)   = [ this_dof-1  this_dof-1  (ELength^2)*Z_val]; % lack of minus sign here is due to E_tot sign in EFIE: - E_scat + E_tot = E_inc 
                Z_rowcolval(Znum+1,:) = [ this_dof  this_dof  0 ];
                %                 
% disp('');
            end
        end
    end




end

 Z_rowcolval(1:1:end,3) = Z_rowcolval(1:1:end,3)*size(interelem_VsrcZload,1); % Normalise by number of vertices