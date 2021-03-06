function [U_Mat] = SelectDOFMBF_FO_2(numVertices, numNodes,numMBF,triangle_blah, U_Mat)

% numDOFS = length(basis_supports);
% U_Mat = zeros(numDOFS,numNodes*3);




phi = 360/numVertices;
% U_vec(triangle_blah(3:
for phi_var = 0:2
    col = 4;
    col_iter = 1;
    total_dofs_selected = 0;
    phi_step = -1;
    iter = -1;
    for i = 1:length(triangle_blah)
        
        if ((i~=(3+(4*iter))) && (i~=4+(4*iter)))
            if (mod(i,2) == 0)
                iter = iter + 1;
                phi_step = phi_step + 2;
            end
            total_dofs_selected = total_dofs_selected + 1;
            
            if total_dofs_selected == (col_iter*(numVertices*2))+1
                
                col_iter = col_iter + 1;
                if (col_iter > numNodes)
                    col = col +numMBF-3;
                else
                    col = col +numMBF;
                end
            end
            
            switch phi_var
                case 0
                    U_Mat(triangle_blah(i,7), col) = 1;
                    U_Mat(triangle_blah(i,8), col) = 1;
                case 1
                    U_Mat(triangle_blah(i,7), col+1) = sind(phi*(phi_step+1));
                    U_Mat(triangle_blah(i,8), col+1) = sind(phi*phi_step);
                case 2
                    U_Mat(triangle_blah(i,7), col+2) = cosd(phi*(phi_step+1));
                    U_Mat(triangle_blah(i,8), col+2) = cosd(phi*phi_step);
            end
        end
        
    end
end