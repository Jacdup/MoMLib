function [U_Mat] = SelectDOFMBF_Circ(basis_supports, numVertices)
% This function selects and numbers the MBF DOFs at each node in the
% contour


numDOFS = length(basis_supports);
% MBF_DOF_Number = 1;
% total_DOFS_selected = 0;
% u_row = 1;
% phi_step = 0;
% iter =1;
U_Vec = zeros(numDOFS,1);
U_Vec_1 = zeros(numDOFS,1);
% U_Vec_sin = zeros(numDOFS,1);
% U_Vec_cos = zeros(numDOFS,1);
% phi_var = 0;
% U_Mat = zeros(numDOFS, numDOFS);
% ----------------------------------------------------------------
% Parse through geometry and assign MBF DOF numbers to each node.
% The row index of the U_Vec is the previous DOF number.
% ----------------------------------------------------------------
phi = 360/numVertices;
for phi_var = 0:2 % All 3 MBFs
    total_DOFS_selected = 0;
    total_DOFS_selected_1 = 0;
    phi_step = -1;
    phi_step_1 = -1;
    MBF_DOF_Number = 1+phi_var;
    MBF_DOF_Number_1 = 4+phi_var;
    iter =1 ;
    iter_1 = 1;
    u_row = 1;
    u_row_1 = 1;
    for i = 1:numDOFS % Total number of DOFS
        if (basis_supports(i,2) - basis_supports(i,1) == numVertices) % If edge sits on a contour node point
            phi_step = phi_step + 1;
            total_DOFS_selected = total_DOFS_selected + 1;
            if total_DOFS_selected == ((iter)*numVertices) % This is the column iterator
                U_Mat(u_row:i,MBF_DOF_Number) = U_Vec(u_row:i);
%                 U_Mat(u_row:i,MBF_DOF_Number+1) = U_Vec_sin(u_row:i);
%                 U_Mat(u_row:i,MBF_DOF_Number+2) = U_Vec_cos(u_row:i);
                MBF_DOF_Number = MBF_DOF_Number + 6;
                iter = iter + 1;
                u_row = i;
%                 U_Vec = zeros(numDOFS,1);
                phi_step = 0;
            end
    %         U_Vec(i) = MBF_DOF_Number;
%             U_Vec(i) = 1; % Change this if want to sample a sinusoid (?)
            switch phi_var
                case 0
                    U_Vec(i) = 1;
                case 1 
                    U_Vec(i) = sind(phi_step*phi);
                case 2
                    U_Vec(i) = cosd(phi_step*phi);
            end
%             U_Vec_sin(i) = sind(phi_step*phi);
%             U_Vec_cos(i) = cosd(phi_step*phi);


        else
           % Add the circumferential MBF in here
            phi_step_1 = phi_step_1 + 1;
            total_DOFS_selected_1 = total_DOFS_selected_1 + 1;
            if total_DOFS_selected_1 == ((iter_1)*numVertices) % This is the column iterator
                U_Mat(u_row_1:i,MBF_DOF_Number_1) = U_Vec_1(u_row_1:i);
                MBF_DOF_Number_1 = MBF_DOF_Number_1 + 6;
                iter_1 = iter_1 + 1;
%                 U_Vec = zeros(numDOFS,1);
                 u_row_1 = i;
                phi_step_1 = 0;
            end
    %         U_Vec(i) = MBF_DOF_Number;
%             U_Vec(i) = 1; % Change this if want to sample a sinusoid (?)
            switch phi_var
                case 0
                    U_Vec_1(i) = 1;
                case 1 
                    U_Vec_1(i) = sind(phi_step_1*phi);
                case 2
                    U_Vec_1(i) = cosd(phi_step_1*phi);
            end
            
        end
    end
    U_Mat(u_row:i,MBF_DOF_Number) = U_Vec(u_row:i);
    U_Mat(u_row_1:i,MBF_DOF_Number_1) = U_Vec_1(u_row_1:i);
end

% U_Mat(u_row:i,MBF_DOF_Number+1) = U_Vec_sin(u_row:i);
% U_Mat(u_row:i,MBF_DOF_Number+2) = U_Vec_cos(u_row:i);






        