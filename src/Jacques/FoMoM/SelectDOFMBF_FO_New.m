function [U_Mat, DOF_mat1, DOF_mat2, DOF_mat3] = SelectDOFMBF_FO_New(mesh_data, dof_data, numVertices ,numMBF, numNodes, triangle_blah, endCap)

% -------------------------------------------------------------------------
% Init
% -------------------------------------------------------------------------
phi = 360/numVertices;
vert_num = (2*numVertices)-1;
if endCap;   numNodes_new = numNodes + 2; endCapExclude = (2*numVertices); % exclude last (2*numVerttices) from i/row assignment
else;    numNodes_new = numNodes; endCapExclude = 0; end
numDofs = size(dof_data.basis_supports,1);
numMBFNodes = numVertices*(numNodes_new);
DOF_mat1 = zeros(numVertices*2,numNodes_new);
DOF_mat2 = zeros(numVertices*2,numNodes_new);
DOF_mat3 = zeros(numVertices*2,numNodes_new);
% X1       = zeros(2,numMBFNodes);
% X2       = zeros(2,numMBFNodes);
% X3       = zeros(2,numMBFNodes);
U_Mat = zeros(numDofs, numNodes_new*numMBF);
% U_Mat = zeros(numDOFS,numNodes*numMBF);

len_tri_mat = (length(triangle_blah)-vert_num-1-endCapExclude);

% -------------------------------------------------------------------------
% Set up MBF matrix
% -------------------------------------------------------------------------
% MBF matrix is the analytical MBF
% MBF_mat has the value of the MBF at each contour node point
% num_nodes x [nodes,constant,sin,cos]
numMBFNodes_new = numVertices*(numNodes+2);
sin_mat = sind(phi*(0:(numMBFNodes_new-1)));
cos_mat = cosd(phi*(0:(numMBFNodes_new-1)));
ones_mat = (ones(numMBFNodes_new,1));
% if endCap == 1
    contour_nodes = (1:numMBFNodes_new)';
%     contour_nodes = (triangle_blah(1,1):(numVertices*(numNodes+1)))'; % All the nodes associated with the analytical MBF
% else
%     contour_nodes = (triangle_blah(2,2):(numVertices*(numNodes_new+1)))'; % All the nodes associated with the analytical MBF
% end
MBF_mat = [contour_nodes,ones_mat,sin_mat',cos_mat'];

% DOF_Mat = ismember(triangle_blah(:,1:3),MBF_mat(:,1)); % Find which triangles contain the node
% DOF_Mat = DOF_Mat .* triangle_blah(:,7:9);
% DOF_Mat = find(triangle_blah(:,1:3) == MBF_mat(:,1));

row = -1;
iter = 1;
col = 1;
% Select DOFs associated with the MBF
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
% -------------------------------------------------------------------------
% -----------------------------Endcap stuff--------------------------------
% -------------------------------------------------------------------------
if endCap == 1
    row1 = row;
    row2 = row;
    for i = 1:2:(vert_num) % triangle_blah is reduced matrix
        row1 = row1+ 2;
        last = 0;
        if (i == vert_num)
            last = 1;
        end

        DOF_mat3(row1:row1+1, col) = [triangle_blah(i,7+last); triangle_blah(i,13+last)];
        DOF_mat1(row1:row1+1, col) = [triangle_blah(i,9); triangle_blah(i,15)];
    end
    for i = (length(triangle_blah)-vert_num-endCapExclude+1):2:length(triangle_blah)-vert_num-1
        row2 = row2 + 2;   
        last = 0;
            if i == length(triangle_blah)-vert_num-1
                last = 1;
            end
            DOF_mat1(row2:row2+1, col+1) = [triangle_blah(i,7); triangle_blah(i,13)];
            DOF_mat2(row2:row2+1, col+1) = [triangle_blah(i,9-last); triangle_blah(i,15-last)];
    end
    DOF_mat1 = [circshift(DOF_mat1(:,1:end-1), [0 1]), DOF_mat1(:,end)]; % Swap columns, to ascending
    DOF_mat2 = [circshift(DOF_mat2(:,1:end-1), [0 1]), DOF_mat2(:,end)];
    DOF_mat3 = [circshift(DOF_mat3(:,1:end-1), [0 1]), DOF_mat3(:,end)];
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

if endCap == 0
    lim = length(edge_vecs_1);
else
    lim = (length(edge_vecs_1) - (numVertices));
end
theta_1      =  abs(90 - acosd(dot(edge_vecs_1(1:lim,:),edge_vecs_2,2)./(vecnorm(edge_vecs_1(1:lim,:),2,2).*vecnorm(edge_vecs_2,2,2))));
theta_2      =  abs(90 - acosd(dot(edge_vecs_1(1:lim,:),edge_vecs_3,2)./(vecnorm(edge_vecs_1(1:lim,:),2,2).*vecnorm(edge_vecs_3,2,2))));

% quiver3(mesh_data.node_coords(edge_nodes_2(:,1),1),mesh_data.node_coords(edge_nodes_2(:,1),2),mesh_data.node_coords(edge_nodes_2(:,1),3),edge_vecs_2(:,1),edge_vecs_2(:,2),edge_vecs_2(:,3))


% Take care of last elements in ring, where local edges are switched
if endCap == 0
%     orientation_vec = (edge_nodes_1 == MBF_mat(edge_nodes_1(1,1): length(edge_nodes_1(:,1)),1))';
    orientation_vec = (edge_nodes_1 == MBF_mat(numVertices+1:end-numVertices,1))'; 
%     orientation_vec = (edge_nodes_1 == MBF_mat(:,1))';
min = 0;
else
    min = numVertices;
    orientation_vec = (edge_nodes_1 == MBF_mat(:,1))';
end
Rho = [1,1;1,-1];
Rho = repmat(Rho,1,1,numMBFNodes);
Rho2 = Rho;
temp = (orientation_vec(2,:) == 1); 
Rho2(:,:,temp(:) == 1) = Rho2(:,:,temp(:) == 1).*[1,-1;1,-1]; % Change minus side when temp == 1
if endCap == 1
    Rho(:,:,1:numVertices) =  Rho(:,:,1:numVertices).*[-1,-1;-1,-1];
else
%      Rho(:,:,1:numVertices) =  Rho(:,:,1:numVertices).*[-1,-1;-1,-1];
%     Rho = Rho2;
%     Rho(:,:,temp(:) == 1) = Rho(:,:,temp(:) == 1).*[1,-1;1,-1]; % Change minus side when temp == 1
end
for MBF_num = 1:3

    B1(1:2,:) = [MBF_mat(edge_nodes_1(:,1),1+MBF_num),MBF_mat(edge_nodes_1(:,2),1+MBF_num)]';
    B2(1:2,:) = [zeros(lim,1),MBF_mat(edge_nodes_2(:,2),1+MBF_num)]';
    B3(1:2,:) = [MBF_mat(edge_nodes_3(:,1),1+MBF_num),zeros(lim,1)]'; 
     B2 = B2(:,:) .* [sind(theta_1)';sind(theta_1)'];
     B3 = B3(:,:) .* [sind(theta_2)';sind(theta_2)'];
     
    for i = 1:numMBFNodes
        X1(:,i) = Rho(:,:,i)\B1(:,i);
    end

    for i = 1:numMBFNodes-min
        X2(:,i) = Rho2(:,:,i)\B2(:,i);
        X3(:,i) = Rho2(:,:,i)\B3(:,i);
    end

    col_iter = 1;
    node2 = 1;
    node3 = 1;
    for MBF_node = 1:numNodes_new
        col_index = col_iter + (MBF_num-1);
        col_iter = col_iter + numMBF;

               U_Mat(DOF_mat1(1:2:end,MBF_node),col_index) = X1(1,(numVertices*(MBF_node-1))+1:(numVertices*MBF_node)); % RWG
               U_Mat(DOF_mat1(2:2:end,MBF_node),col_index) = X1(2,(numVertices*(MBF_node-1))+1:(numVertices*MBF_node)); % Linear
               
               % Skip the zero columns
               if DOF_mat2(1,MBF_node) ~= 0
                   xdom = (numVertices*(node2-1))+1:(numVertices*node2);
                   U_Mat(DOF_mat2(1:2:end,MBF_node),col_index) = X2(1,xdom); % RWG
                   U_Mat(DOF_mat2(2:2:end,MBF_node),col_index) = X2(2,xdom); % Linear
                   % Move to next X domain
                   node2 = node2 + 1;
               end
               if DOF_mat3(1,MBF_node) ~= 0
                    xdom = (numVertices*(node3-1))+1:(numVertices*node3);
                    U_Mat(DOF_mat3(1:2:end,MBF_node),col_index) = X3(1,xdom); % RWG
                    U_Mat(DOF_mat3(2:2:end,MBF_node),col_index) = X3(2,xdom); % Linear
                    node3 = node3 + 1;
               end

    end
end




% -----------------------------Endcap stuff--------------------------------
% Add the first endcap MBF
if endCap == 2 % temporary to not loop this
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