function [U_Mat] = SelectDOFMBF_FO_New(mesh_data, dof_data, numVertices ,numMBF,triangle_blah, endCap)

% numDOFS = length(basis_supports);
% triangle_blah = mesh_data.reduced
% U_Mat = zeros(numDOFS,numNodes*numMBF);
dof_data.dofs_to_edges = dof_data.dofs_to_edges(1:2:end,:); % Every 2nd element
if endCap == 1
    
    % exclude last (2*numVerttices) from i/row assignment
    endCapExclude = (2*numVertices);
    
else
    endCapExclude = 0;
end

total_dofs_selected = 0; %TODO: divide this by 3 to get one column DOFs
phi = 360/numVertices;
vert_num = (2*numVertices)-1;
% U_vec(triangle_blah(3:
sin_mat = sind(phi*(1:numVertices));
cos_mat = cosd(phi*(1:numVertices));
tri_num = -3;

for phi_var = 0:2
    col = 1;
    %     col_iter = 1;
    
    phi_step = -1;
    iter = 1;
    
    for i = 1:2:(length(triangle_blah)-vert_num-1-endCapExclude) % triangle_blah is reduced matrix
%         tri_num = tri_num + 4; % Every reduced matrix increment by 2 means the original is incremented by 4
        phi_step = phi_step + 1;
        total_dofs_selected = total_dofs_selected + 1;
        % Assign the three edges (RWG and FO) of each azimuth point
        % This is specific to how the DOF's are assigned in the
        % preprocessing step
        if i == ((vert_num+1)*iter)-1
            last = 1;
            iter = iter + 1;
        else
            last = 0;
        end
        
%         if tri_num > 192
%             test = 1;
%         end
         edge = 1 + last;% Local edge
        Rho = [1,-1;1,1];

        switch phi_var
            case 0
                U_Mat_num = 1;
                X(1) = 1;
                X(2) = 0; % For constant, set linear to 0.
                col_iter =  0;
            case 1
%                 U_Mat_num = sind(phi*phi_step);
                col_iter =  1;
                
                trig_mat = sin_mat;
%                 B = [sin_mat(mod(phi_step, numVertices)+1); sin_mat(mod(phi_step+1,numVertices)+1)];
%                 X = Rho\B; % Calculate coefficients for both RWG and Linear 
            case 2
%                 U_Mat_num = cosd(phi*phi_step);
                trig_mat = cos_mat;
%                 B = [cos_mat(mod(phi_step, numVertices)+1); cos_mat(mod(phi_step+1,numVertices)+1)];
%                 X = Rho\B; % Calculate coefficients for both RWG and Linear 
                col_iter =  2;
        end

       
        % Node edge
        if phi_var == 0
            X(1) = 1;
            X(2) = 0 ;
        else
%             [Rho] = CalculateBasisFunction(mesh_data, triangle_blah, 1, i+1);
             B = [trig_mat(mod(phi_step, numVertices)+1); trig_mat(mod(phi_step+1,numVertices)+1)];
             X = Rho\B;
        end
        U_Mat(triangle_blah(i+1,7),col+col_iter) = X(1);
        U_Mat(triangle_blah(i+1,13),col+col_iter) = X(2);
        
        % First diagonal edge
        % TODO: check which is lower node (for 0)
        if phi_var == 0
%                 if (triangle_blah(i+last,7+last) == 5) || (triangle_blah(i+last,7+last) == 6)
%                  test = 1;
%                 end
%              [Rho] = CalculateBasisFunction(mesh_data, triangle_blah, edge, i+last);
               

%              edge_nodes =
%              mesh_data.edges(dof_data.dofs_to_edges(i+last),:); % This
%              gives nodes associated with edge. edge_nodes(1) gives lower
%              node.
             B = [0; 1];
             X = Rho\B;
%             X(1) = 1;
%             X(2) = 0 ;
        else
%             [Rho] = CalculateBasisFunction(mesh_data, triangle_blah, edge, i+last);
            B = [0;trig_mat(mod(phi_step, numVertices)+1)];  
            X = Rho\B;
        end
        U_Mat(triangle_blah(i+last,7+last),col+col_iter) = X(1);
        U_Mat(triangle_blah(i+last,13+last), col+col_iter) = X(2);
        
        % Second diagonal edge
        if phi_var == 0
%             
%             [Rho] = CalculateBasisFunction(mesh_data, triangle_blah, edge, i+vert_num+last);
            
             B = [1; 0];
             X = Rho\B;
%             X(1) = 1;
%             X(2) = 0 ;
        else
%             [Rho] = CalculateBasisFunction(mesh_data, triangle_blah, edge, i+vert_num+1);
            B = [trig_mat(mod(phi_step, numVertices)+1);0 ];
            X = Rho\B;
        end
        U_Mat(triangle_blah(i+vert_num+1,7+last),col+col_iter) = X(1);
        U_Mat(triangle_blah(i+vert_num+1,13+last),col+col_iter) = X(2);
        
        
        if last == 1
            col = col + numMBF;
        end
        
        
    end
    U_Mat(max(max(triangle_blah))-1,:) = 0;
    U_Mat(max(max(triangle_blah)),:) = 0;
end

% Add the first endcap MBF
if endCap == 1
for phi_var = 0:2
    phi_step = -1;
    phi_step1 = -1;

    
    
    
    for i = 1:2:(vert_num) % triangle_blah is reduced matrix
        last = 0;
        phi_step = phi_step + 1;
        if (i == vert_num)
            last = 1;
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
        
        U_Mat(triangle_blah(i,7+last),col+col_iter) = U_Mat_num;
        U_Mat(triangle_blah(i,13+last),col+col_iter) = U_Mat_num;
        U_Mat(triangle_blah(i,9),col+col_iter) = U_Mat_num;
        U_Mat(triangle_blah(i,15),col+col_iter) = U_Mat_num;
    end
    
    % 2nd endcap MBF
    for i = (length(triangle_blah)-vert_num-endCapExclude+1):2:length(triangle_blah)-vert_num-1
        last = 0;
        phi_step1 = phi_step1 + 1;
        if i == length(triangle_blah)-vert_num-1
            last = 1;
            
        end
        switch phi_var
            case 0
                U_Mat_num = 1;
                col_iter =  0;
            case 1
                U_Mat_num = sind(phi*phi_step1);
                col_iter =  1;
            case 2
                U_Mat_num = cosd(phi*phi_step1);
                col_iter =  2;
        end
        
        U_Mat(triangle_blah(i,7),col+col_iter+numMBF) = U_Mat_num;
        U_Mat(triangle_blah(i,13),col+col_iter+numMBF) = U_Mat_num;
        U_Mat(triangle_blah(i,9-last),col+col_iter+numMBF) = U_Mat_num;
        U_Mat(triangle_blah(i,15-last),col+col_iter+numMBF) = U_Mat_num;
    end
    
    
end



end