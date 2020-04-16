function [U_Mat] = SelectDOFMBF_FO_New(mesh_data, dof_data, numVertices ,numMBF, numNodes, triangle_blah, endCap)


total_dofs_selected = 0; %TODO: divide this by 3 to get one column DOFs
phi = 360/numVertices;
vert_num = (2*numVertices)-1;
numDofs = size(dof_data.basis_supports,1);
numMBFNodes = numVertices*numNodes;
DOF_mat1 = zeros(numVertices*2,numNodes);
DOF_mat2 = zeros(numVertices*2,numNodes);
DOF_mat3 = zeros(numVertices*2,numNodes);
X1       = zeros(2,numMBFNodes);
X2       = zeros(2,numMBFNodes);
X3       = zeros(2,numMBFNodes);
% triangle_blah = mesh_data.reduced
% U_Mat = zeros(numDOFS,numNodes*numMBF);

if endCap == 1
    % exclude last (2*numVerttices) from i/row assignment
    endCapExclude = (2*numVertices);
else
    endCapExclude = 0;
end
len_tri_mat = (length(triangle_blah)-vert_num-1-endCapExclude);

% Set up MBF matrix
% MBF matrix is the analytical MBF
% MBF_mat has the value of the MBF at each contour node point
% num_nodes x [nodes,constant,sin,cos]
sin_mat = sind(phi*(1:(numMBFNodes)));
cos_mat = cosd(phi*(1:(numMBFNodes)));
contour_nodes = (triangle_blah(2,2):(numVertices*(numNodes+1)))'; % All the nodes associated with the analytical MBF
MBF_mat = [contour_nodes,ones(numMBFNodes, 1),sin_mat',cos_mat'];

% DOF_Mat = ismember(triangle_blah(:,1:3),MBF_mat(:,1)); % Find which triangles contain the node
% DOF_Mat = DOF_Mat .* triangle_blah(:,7:9);
% DOF_Mat = find(triangle_blah(:,1:3) == MBF_mat(:,1));

row = -1;
iter = 1;
col = 1;
for i = 1:2:len_tri_mat
    row = row + 2;
    if i == ((vert_num+1)*iter)-1
            last = 1;
            iter = iter + 1;
    else
            last = 0;
    end
    
    % Create matrix containing DOFs of the edges
    DOF_mat1(row:row+1,col) = [triangle_blah(i+1,7);triangle_blah(i+1,13)];
    DOF_mat2(row:row+1,col) = [triangle_blah(i+last,7+last);triangle_blah(i+last,13+last)];
    DOF_mat3(row:row+1,col) = [triangle_blah(i+vert_num+1,7+last);triangle_blah(i+vert_num+1,13+last)];
    
    if last == 1
        row = -1;
        col = col + 1;
    end
    
end
temp1        = nonzeros(DOF_mat1); % Create temporary column vector
temp2        = nonzeros(DOF_mat2);
temp3        = nonzeros(DOF_mat3);
edge_nodes_1 = mesh_data.edges(dof_data.dofs_to_edges(temp1(1:2:end,1)),:); % Edges on contour
edge_nodes_2 = mesh_data.edges(dof_data.dofs_to_edges(temp2(1:2:end,1)),:); % First diagonal
edge_nodes_3 = mesh_data.edges(dof_data.dofs_to_edges(temp3(1:2:end,1)),:); % Second diagonal

edge_vecs_1  = mesh_data.node_coords(edge_nodes_1(:,1),:)-mesh_data.node_coords(edge_nodes_1(:,2),:);
edge_vecs_2  = mesh_data.node_coords(edge_nodes_2(:,1),:)-mesh_data.node_coords(edge_nodes_2(:,2),:);
edge_vecs_3  = mesh_data.node_coords(edge_nodes_3(:,1),:)-mesh_data.node_coords(edge_nodes_3(:,2),:);

theta_1      =  abs(90 - acosd(dot(edge_vecs_1,edge_vecs_2,2)./(vecnorm(edge_vecs_1,2,2).*vecnorm(edge_vecs_2,2,2))));
theta_2      =  abs(90 - acosd(dot(edge_vecs_1,edge_vecs_3,2)./(vecnorm(edge_vecs_1,2,2).*vecnorm(edge_vecs_3,2,2))));

% quiver3(mesh_data.node_coords(edge_nodes_2(:,1),1),mesh_data.node_coords(edge_nodes_2(:,1),2),mesh_data.node_coords(edge_nodes_2(:,1),3),edge_vecs_2(:,1),edge_vecs_2(:,2),edge_vecs_2(:,3))
% theta_3      = (pi/2) - acos(dot(edge_vecs_1,edge_vecs_2,2)/(norm(edge_vecs_1)*norm(edge_vecs_2)));

B1 = (edge_nodes_1 == MBF_mat(:,1))'; % TODO, this should be all ones, since edge_nodes_1 contain all nodes in MBF_mat
B2 = (edge_nodes_2 == MBF_mat(:,1))'; 
B3(1:numMBFNodes,1:2) = [ones(numMBFNodes,1), zeros(numMBFNodes,1)];% TODO, do this generically. I can do this now because I know a priori how
% this should look
B3 = logical(B3');




Rho = [1,1;1,-1];
Rho = repmat(Rho,1,1,numMBFNodes);
temp = (B1(2,:) == 1); 
Rho(:,:,temp(:) == 1) = Rho(:,:,temp(:) == 1).*[1,-1;1,-1]; % Change minus side when temp == 1

B1 = B1(:,:) .* [MBF_mat(:,2),MBF_mat(:,2)]';
B2 = B2(:,:) .* [MBF_mat(:,2).*sind(theta_1),MBF_mat(:,2).*sind(theta_1)]';
B3 = B3(:,:) .* [MBF_mat(:,2).*sind(theta_2),MBF_mat(:,2).*sind(theta_2)]';

for i = 1:numMBFNodes
    X1(:,i) = Rho(:,:,i)\B1(:,i);
    X2(:,i) = Rho(:,:,i)\B2(:,i);
    X3(:,i) = Rho(:,:,i)\B3(:,i);
end


% DOF_mat(1:2:end,:) is the RWG
% DOF_mat(2:2:end,:) is linear
U_Mat = zeros(numDofs, numNodes);
for MBF_node = 1:numNodes
    
%     X1 = pagefun(@mldivide,gpuArray(Rho),gpuArray(B1));
%     U_Mat(DOF_mat1(1:2:end,MBF_node),MBF_node) = X1(1,(numVertices*(MBF_node-1))+1:(numVertices*MBF_node)); % RWG 
    U_Mat(DOF_mat1(1:2:end,MBF_node),MBF_node) = 1;
    U_Mat(DOF_mat1(2:2:end,MBF_node),MBF_node) = 0; % Linear

    U_Mat(DOF_mat2(1:2:end,MBF_node),MBF_node) = X2(1,(numVertices*(MBF_node-1))+1:(numVertices*MBF_node)); % RWG 
    U_Mat(DOF_mat2(2:2:end,MBF_node),MBF_node) = X2(2,(numVertices*(MBF_node-1))+1:(numVertices*MBF_node)); % Linear

    U_Mat(DOF_mat3(1:2:end,MBF_node),MBF_node) = X3(1,(numVertices*(MBF_node-1))+1:(numVertices*MBF_node)); % RWG 
    U_Mat(DOF_mat3(2:2:end,MBF_node),MBF_node) = X3(2,(numVertices*(MBF_node-1))+1:(numVertices*MBF_node)); % Linear
    
 
end

% U_Mat(DOF_mat1(numVertices+1:numVertices*2,2)) = MBF_mat(:,2); % Assign the vals. TODO, multiply by B
% U_Mat(DOF_mat2(numVertices+1:numVertices*2,2)) = MBF_mat(:,2);
% U_Mat(DOF_mat3(numVertices+1:numVertices*2,2)) = MBF_mat(:,2);


% dof_data.dofs_to_edges = dof_data.dofs_to_edges(1:2:end,:); % Every 2nd element

% for phi_var = 0:2
%     col = 1;
%     %     col_iter = 1;
%     
%     phi_step = -1;
%     iter = 1;
%     row = 0;
%     for i = 1:2:len_tri_mat % triangle_blah is reduced matrix
%         %         tri_num = tri_num + 4; % Every reduced matrix increment by 2 means the original is incremented by 4
%         phi_step = phi_step + 1;
%         total_dofs_selected = total_dofs_selected + 1;
%         row = row + 1;
%         % Assign the three edges (RWG and FO) of each azimuth point
%         % This is specific to how the DOF's are assigned in the
%         % preprocessing step
%         if i == ((vert_num+1)*iter)-1
%             last = 1;
%             iter = iter + 1;
%         else
%             last = 0;
%         end
%         
%         %         if tri_num > 192
%         %             test = 1;
%         %         end
%         edge = 1 + last;% Local edge
%         
%         
%         switch phi_var
%             case 0
%                 U_Mat_num = 1;
%                 %                 X(1) = 1;
%                 %                 X(2) = 0; % For constant, set linear to 0.
%                 col_iter =  0;
%             case 1
%                 %                 U_Mat_num = sind(phi*phi_step);
%                 col_iter =  1;
%                 
%                 %                 trig_mat = sin_mat;
%                 %                 U_Mat_num = trig_mat(mod(phi_step+1, numVertices)+1);
%                 %                 B = [sin_mat(mod(phi_step, numVertices)+1); sin_mat(mod(phi_step+1,numVertices)+1)];
%                 %                 X = Rho\B; % Calculate coefficients for both RWG and Linear
%             case 2
%                 %                 U_Mat_num = cosd(phi*phi_step);
%                 %                 trig_mat = cos_mat;
%                 %                 U_Mat_num = trig_mat(mod(phi_step+1, numVertices)+1);
%                 %                 B = [cos_mat(mod(phi_step, numVertices)+1); cos_mat(mod(phi_step+1,numVertices)+1)];
%                 %                 X = Rho\B; % Calculate coefficients for both RWG and Linear
%                 col_iter =  2;
%         end
%         Rho          = [1,1;1,-1];
%         edge_nodes_1 = mesh_data.edges(dof_data.dofs_to_edges(triangle_blah(i+1,7)),:); % This
%         %              gives nodes associated with edge. edge_nodes(1) gives lower
%         %              node.
%         edge_nodes_2 = mesh_data.edges(dof_data.dofs_to_edges(triangle_blah(i+last,7+last)),:);
%         edge_nodes_3 = mesh_data.edges(dof_data.dofs_to_edges(triangle_blah(i+vert_num+1, 7+last)),:);
%         if edge_nodes_1(2) == MBF_mat(row, 1) % If the higher number is a contour node
%             % TODO: check this
%             Rho = [1,-1;1,1];
%         end
%         
%         if phi_var == 0
%             X1(1) = 1;
%             X1(2) = 0; % Linear is zero for constant case
%         else
%             next_val = mod(row+1, numVertices) ;
%             if next_val == 0
%                 next_val = next_val + 1;
%             end
%             if edge_nodes_1(1) == MBF_mat(row,1)
%                 % TODO: check this
%                 B1 = [MBF_mat(row,2+phi_var); MBF_mat(next_val,2+phi_var)];
%             else
%                 B1 = [MBF_mat(mod(row+1,numVertices)+1,2+phi_var);MBF_mat(row,2+phi_var) ];
%             end
%             X1 = Rho\B1;
%         end
%         
%         if edge_nodes_2(1) == MBF_mat(row,1)
%             B2 = [MBF_mat(row,2+phi_var);0];
%         else %elseif edge_nodes_2(2) = MBF_mat(row,1) % This should also be true
%             B2 = [0;MBF_mat(row,2+phi_var)];
%         end
%         
%         if edge_nodes_3(1) == MBF_mat(row,1)
%             B3 = [MBF_mat(row,2+phi_var);0];
%         else %elseif edge_nodes_3(2) = MBF_mat(row,1) % This should also be true
%             B3 = [0;MBF_mat(row,2+phi_var)];
%         end
%         
%         
%         % Node edge
%         U_Mat(triangle_blah(i+1,7),col+col_iter) = X1(1);
%         U_Mat(triangle_blah(i+1,13),col+col_iter) = X1(2);
%         
%         % First diagonal edge
%         %          B = [0;trig_mat(mod(phi_step+1, numVertices)+1)];
%         X2 = Rho\B2;
%         U_Mat(triangle_blah(i+last,7+last),col+col_iter) = X2(1);
%         U_Mat(triangle_blah(i+last,13+last), col+col_iter) = X2(2);
%         
%         % Second diagonal edge
%         %         B = [trig_mat(mod(phi_step, numVertices)+1);0 ];
%         X3 = Rho\B3;
%         U_Mat(triangle_blah(i+vert_num+1,7+last),col+col_iter) = X3(1);
%         U_Mat(triangle_blah(i+vert_num+1,13+last),col+col_iter) = X3(2);
%         
%         
%         if last == 1
%             col = col + numMBF;
%         end
%         
%         
%     end
%     U_Mat(max(max(triangle_blah))-1,:) = 0;
%     U_Mat(max(max(triangle_blah)),:) = 0;
% end

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