function [U_Mat] = SelectDOFMBF_2(basis_supports, numVertices, U_Mat, numNodes,numMBF)
% This function selects and numbers the MBF DOFs at each node in the
% contour


numDOFS = length(basis_supports);
U_Vec = zeros(numDOFS,1);

% U_Mat = zeros(numDOFS, numDOFS);
% ----------------------------------------------------------------
% Parse through geometry and assign MBF DOF numbers to each node.
% The row index of the U_Vec is the previous DOF number.
% ----------------------------------------------------------------
phi = 360/numVertices;
for phi_var = 0:2 % All 3 MBFs
    total_DOFS_selected = 0;
    phi_step = -1;
    MBF_DOF_Number = (numMBF-2)+phi_var;
    iter =1 ;
    u_row = 1;
    for i = 1:numDOFS % Total number of DOFS
        
        if (basis_supports(i,2) - basis_supports(i,1) ~= numVertices) % If edge sits on a contour node point
            phi_step = phi_step + 1;
            total_DOFS_selected = total_DOFS_selected + 1;
            
            if total_DOFS_selected == ((iter)*numVertices)+1 % This is the column iterator
                U_Mat(u_row:i-1,MBF_DOF_Number) = U_Vec(u_row:i-1);
                
                iter = iter +1;
                if (iter > numNodes)
                    MBF_DOF_Number = MBF_DOF_Number + numMBF-3;
                else
                    MBF_DOF_Number = MBF_DOF_Number + numMBF;
                end
                
                u_row = i;
                phi_step = 0;
            end
            
            switch phi_var
                case 0
                    U_Vec(i) = 1;
                case 1
                    %                     U_Vec(i) = 1;
                    U_Vec(i) = 1*sind(phi_step*phi);
                case 2
                    %                     U_Vec(i) = 1;
                    U_Vec(i) = 1*cosd(phi_step*phi);
            end
        end
    end
    U_Mat(u_row:i,MBF_DOF_Number) = U_Vec(u_row:i);
end

% U_Mat(u_row:i,MBF_DOF_Number+1) = U_Vec_sin(u_row:i);
% U_Mat(u_row:i,MBF_DOF_Number+2) = U_Vec_cos(u_row:i);






