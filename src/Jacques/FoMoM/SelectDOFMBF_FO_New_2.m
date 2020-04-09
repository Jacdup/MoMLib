function [U_Mat] = SelectDOFMBF_FO_New_2(numVertices, numNodes,numMBF,triangle_blah, U_Mat)

% numDOFS = length(basis_supports);
% U_Mat = zeros(numDOFS,numNodes*3);

phi = 360/numVertices;

for phi_var = 0:2
    col = numMBF-2;
    col_iter = 1;
    total_dofs_selected = 0;
    phi_step = -1;
    i_iter = 0;
    for i = 1:2:length(triangle_blah)% Every odd row
        i_iter = i_iter + 1;
        phi_step = phi_step + 2;
        total_dofs_selected = total_dofs_selected + 1;
        
        switch phi_var
            case 0
                U_Mat_num       = 1;
                U_Mat_num_next  = 1;
                col_iter        = 0;
            case 1
                U_Mat_num       = sind(phi*phi_step);
                U_Mat_num_next  = sind(phi*(phi_step+1));
                col_iter        = 1;
            case 2
                U_Mat_num       = cosd(phi*phi_step);
                U_Mat_num_next  = cosd(phi*(phi_step+1));
                col_iter        = 2;
        end
        
        U_Mat(triangle_blah(i,8),col+col_iter)  = U_Mat_num;
        U_Mat(triangle_blah(i,14),col+col_iter) = U_Mat_num;
        U_Mat(triangle_blah(i,7),col+col_iter)  = U_Mat_num_next; % Next angle here
        U_Mat(triangle_blah(i,13),col+col_iter) = U_Mat_num_next;
        
        if i_iter == numVertices % Round we go
            i_iter = 0;
            
            U_Mat(triangle_blah(i,7),col+col_iter)  = U_Mat_num;
            U_Mat(triangle_blah(i,13),col+col_iter) = U_Mat_num;  
            U_Mat(triangle_blah(i,8),col+col_iter)  = U_Mat_num_next; % Next angle here
            U_Mat(triangle_blah(i,14),col+col_iter) = U_Mat_num_next;
            phi_step = -1;
            if (ceil(col/numMBF)) >= numNodes
                col = col + numMBF - 3;
            else
                col = col + numMBF;
            end
        end
        
    end
end