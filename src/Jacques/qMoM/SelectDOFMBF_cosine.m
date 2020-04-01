function [U_Mat] = SelectDOFMBF(basis_supports, numVertices)%filename recieved is a string
% This function selects and numbers the MBF DOFs at each node in the
% contour


numDOFS = length(basis_supports);
MBF_DOF_Number = 1;
total_DOFS_selected = 0;
u_row = 1;
phi_step = 0;

U_Vec = zeros(numDOFS,1);
% U_Mat = zeros(numDOFS, numDOFS);
% ----------------------------------------------------------------
% Parse through geometry and assign MBF DOF numbers to each node.
% The row index of the U_Vec is the previous DOF number.
% ----------------------------------------------------------------
phi = 360/numVertices;
for i = 1:numDOFS % Total number of DOFS
    if (basis_supports(i,2) - basis_supports(i,1) == numVertices)
        phi_step = phi_step + 1;
        total_DOFS_selected = total_DOFS_selected + 1;
        if total_DOFS_selected == (MBF_DOF_Number*numVertices)+1
            U_Mat(u_row:i,MBF_DOF_Number) = U_Vec(u_row:i);
            MBF_DOF_Number = MBF_DOF_Number + 1;
            u_row = i;
            phi_step = 0;
        end
%         U_Vec(i) = MBF_DOF_Number;
%         U_Vec(i) = 1; % Change this if want to sample a sinusoid (?)
        U_Vec(i) = 1;
        
    end
end
U_Mat(u_row:i,MBF_DOF_Number) = U_Vec(u_row:i);






        