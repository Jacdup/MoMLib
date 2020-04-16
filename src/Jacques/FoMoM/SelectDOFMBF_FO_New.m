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
U_Mat = zeros(numDofs, numNodes*numMBF);
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
B1 = ones(2,numMBFNodes);

for MBF_num = 1:numMBF
    
    B1 = B1(:,:) .* [MBF_mat(:,1+MBF_num),circshift(MBF_mat(:,1+MBF_num),-1)]';
    B2 = B2(:,:) .* [MBF_mat(:,1+MBF_num).*sind(theta_1),circshift(MBF_mat(:,1+MBF_num).*sind(theta_1),-1)]';
    B3 = B3(:,:) .* [MBF_mat(:,1+MBF_num).*sind(theta_2),circshift(MBF_mat(:,1+MBF_num).*sind(theta_2),-1)]';
    
    for i = 1:numMBFNodes
        X1(:,i) = Rho(:,:,i)\B1(:,i);
        X2(:,i) = Rho(:,:,i)\B2(:,i);
        X3(:,i) = Rho(:,:,i)\B3(:,i);
    end

    col_iter = 1;
    for MBF_node = 1:numNodes
        col_index = col_iter + (MBF_num-1);
        col_iter = col_iter + numMBF;

        U_Mat(DOF_mat1(1:2:end,MBF_node),col_index) = X1(1,(numVertices*(MBF_node-1))+1:(numVertices*MBF_node)); % RWG
        U_Mat(DOF_mat1(2:2:end,MBF_node),col_index) = X1(2,(numVertices*(MBF_node-1))+1:(numVertices*MBF_node)); % Linear

        U_Mat(DOF_mat2(1:2:end,MBF_node),col_index) = X2(1,(numVertices*(MBF_node-1))+1:(numVertices*MBF_node)); % RWG
        U_Mat(DOF_mat2(2:2:end,MBF_node),col_index) = X2(2,(numVertices*(MBF_node-1))+1:(numVertices*MBF_node)); % Linear
        
        U_Mat(DOF_mat3(1:2:end,MBF_node),col_index) = X3(1,(numVertices*(MBF_node-1))+1:(numVertices*MBF_node)); % RWG
        U_Mat(DOF_mat3(2:2:end,MBF_node),col_index) = X3(2,(numVertices*(MBF_node-1))+1:(numVertices*MBF_node)); % Linear
        
    end
end




% -----------------------------Endcap stuff--------------------------------
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