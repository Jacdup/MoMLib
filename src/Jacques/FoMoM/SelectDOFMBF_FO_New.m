function [U_Mat] = SelectDOFMBF_FO_New(numDOFS, numVertices, numNodes,numMBF,triangle_blah)

% numDOFS = length(basis_supports);

% U_Mat = zeros(numDOFS,numNodes*numMBF);
total_dofs_selected = 0; %TODO: divide this by 3 to get one column DOFs
phi = 360/numVertices;
vert_num = (2*numVertices)-1;
% U_vec(triangle_blah(3:
for phi_var = 0:2
    col = 1;
%     col_iter = 1;
    
    phi_step = -1;
    iter = 1;
    
    for i = 1:2:(length(triangle_blah)-vert_num-1) % triangle_blah is reduced matrix
        
        phi_step = phi_step + 1;
        total_dofs_selected = total_dofs_selected + 1;
        % Assign the three edges (RWG and FO) of each azimuth point
        % TODO: add sine and cosine
        % This is specific to how the DOF's are assigned in the
        % preprocessing step
        if i == ((vert_num+1)*iter)-1
            
            last = 1;
            iter = iter + 1;
            
        else
            last = 0;
        end
        
        switch phi_var
            case 0
                U_Mat_num = 1;
                col_iter =  0;
            case 1
                U_Mat_num = sind(phi*phi_step);
                col_iter =  1;
            case 2
                U_Mat_num = cosd(phi*phi_step);
                col_iter =  2;
        end
%         i
%         if i == 37
%             test = 1;
%         end
        U_Mat(triangle_blah(i+last,7+last),col+col_iter) = U_Mat_num;
        U_Mat(triangle_blah(i+last,13+last), col+col_iter) = U_Mat_num;
        U_Mat(triangle_blah(i+1,7),col+col_iter) = U_Mat_num;
        U_Mat(triangle_blah(i+1,13),col+col_iter) = U_Mat_num;
        U_Mat(triangle_blah(i+vert_num+1,7+last),col+col_iter) = U_Mat_num;
        U_Mat(triangle_blah(i+vert_num+1,13+last),col+col_iter) = U_Mat_num;
        
        if last == 1
            col = col + 3;
        end
        
        
    end
    U_Mat(max(max(triangle_blah))-1,:) = 0;
    U_Mat(max(max(triangle_blah)),:) = 0;
end