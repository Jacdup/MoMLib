function [Z_lump, V_lump] = FinalZmatVvecLumped(Z_rowcolval, V_rowcolval, global_to_redu_dofs, num_obs, num_src, observer_map, source_map)
% This is a simple routime which takes the lumped element (voltage sources
% and lumped impedances) contributions to the global system of linear
% equations and reduces it to the corresponding contributions to the
% reduced system. I.e. the indices are updated and contributions associated
% with dropped rows/cols are dropped.
%
% The end results are stored in sparse matrices which can simply be added
% to the final MoM Z_mat and V_vec.
%
% 2019-12-16: Created. MMB.

% Init:
num_lumped  = size(Z_rowcolval,1); % max (max possible number of lumped V's, max possible number of lumped Z's)
Z_lump      = spalloc(num_obs, num_src, min(num_lumped,num_obs*num_src));
V_lump      = spalloc(num_obs,       1, min(num_lumped,num_obs        ));


% Add lumped circuit element contributions to the reduced Z and V:
for ii = 1:size(Z_rowcolval,1)
    if min(Z_rowcolval(ii,1:2),2) > 0 % this is a valid entry in the global matrix (can be invalid, since Zlump_rowcolval can have rows of unspecified Z values, with indices (-1,-1)
        new_idx = global_to_redu_dofs(1,Z_rowcolval(ii,1:2));
        new_idx = [ observer_map(new_idx(1,1))  source_map(new_idx(1,2)) ];
        if min(new_idx,2) > 0 % this is a valid entry in the reduced matrix
            Z_lump(new_idx(1,1),new_idx(1,2)) = Z_lump(new_idx(1,1),new_idx(1,2)) + Z_rowcolval(ii,3);
        end
    end
end
for ii = 1:size(V_rowcolval,1)
    if V_rowcolval(ii,1) > 0 % this is a valid entry in the global vector (only row index is of interest now --- observers)
        new_idx = global_to_redu_dofs(1,V_rowcolval(ii,1));
        new_idx = observer_map(new_idx);
        if new_idx > 0 % this is a valid entry in the reduced vector
            V_lump(new_idx,1) = V_lump(new_idx,1) + V_rowcolval(ii,3);
        end
    end
end
