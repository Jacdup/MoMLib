function [U_Mat, DOF_mat1, DOF_mat2, DOF_mat3] = MBF_Axial(mesh_data, dof_data, numVertices ,numMBF, numNodes, triangle_blah, connection, cyl_def)


oneEndcap = (cyl_def.firstNode == "endCap" && cyl_def.lastNode ~= "endCap") || (cyl_def.firstNode ~= "endCap" && cyl_def.lastNode == "endCap"); % Only one is an endcap
twoEndcaps = cyl_def.firstNode == "endCap" && cyl_def.lastNode == "endCap";
% -------------------------------------------------------------------------
% Init
% -------------------------------------------------------------------------
phi = 360/numVertices;
vert_num = (2*numVertices)-1;
if twoEndcaps
% if endCap || connection
     % In the connection case, the setup is exactly the same as the endcap,
    % since it is a special case of the endcap
    numNodes_new = numNodes + 2;
    endCapExclude = (2*numVertices); % exclude last (2*numVertices) from i/row assignment
elseif oneEndcap
    numNodes_new = numNodes + 1; 
    endCapExclude = numVertices;
else
    numNodes_new = numNodes; 
    endCapExclude = 0;
end

numDofs = size(dof_data.basis_supports,1);
numMBFNodes = numVertices*(numNodes_new);
DOF_mat1 = zeros(numVertices*2,numNodes_new);
DOF_mat2 = zeros(numVertices*2,numNodes_new);
DOF_mat3 = zeros(numVertices*2,numNodes_new);
X1       = zeros(2,numMBFNodes);
X2       = zeros(2,numMBFNodes);
X3       = zeros(2,numMBFNodes);
U_Mat = zeros(numDofs, numNodes_new*numMBF);

len_tri_mat = (length(triangle_blah)-vert_num-1-endCapExclude);

Rho = [1,1;1,-1]; % This is the matrix of RWG and linear components at the two edge nodes
Rho = repmat(Rho,1,1,numMBFNodes);
Rho2 = Rho;

% -------------------------------------------------------------------------
% Set up MBF matrix
% -------------------------------------------------------------------------
% MBF matrix is the analytical MBF
% MBF_mat has the value of the MBF at each contour node point
% num_nodes x [nodes,constant,sin,cos]
% numMBFNodes_new = numVertices*(numNodes+4); % For coupling
numMBFNodes_new = numVertices*(numNodes+2); % This is only for filling the MBF_mat
sin_mat = sind(phi*(0:(numMBFNodes_new-1)));
cos_mat = cosd(phi*(0:(numMBFNodes_new-1)));
ones_mat = (ones(numMBFNodes_new,1));
contour_nodes = (1:numMBFNodes_new)';

MBF_mat = [contour_nodes,ones_mat,sin_mat',cos_mat'];
% MBF_mat = MBF_mat(1:end-36,:); % Temporary for coupling

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
for i = numVertices:numVertices:length(Rho)-1 % DOFs at the end of each node have switched signs
    Rho2(:,:,i) = Rho2(:,:,i) .* [1,-1;1,-1];
end

% Temporary for coupling
% DOF_mat1(:,32) = [];
% DOF_mat2(:,32) = [];
% DOF_mat3(:,32) = [];
% -------------------------------------------------------------------------
% -----------------------------Endcap stuff--------------------------------
% -------------------------------------------------------------------------
if oneEndcap || twoEndcaps %|| connection 
    row1 = -1; % Fill up a new column, from row 1
    row2 = -1;
    
    if cyl_def.firstNode == "endCap"
        for i = 1:2:(vert_num) % First endcap, or connection triangles
            row1 = row1+ 2;
            last = 0;
            if (i == vert_num)
                last = 1;
            end

            DOF_mat3(row1:row1+1, col) = [triangle_blah(i,7+last); triangle_blah(i,13+last)];
            DOF_mat1(row1:row1+1, col) = [triangle_blah(i,9); triangle_blah(i,15)];
        end
    end
    if cyl_def.lastNode == "endCap"
        for i = (length(triangle_blah)-vert_num-endCapExclude+1):2:length(triangle_blah)-endCapExclude % Second endcap, or connection triangles
            row2 = row2 + 2;
            last = 0;
            if i == length(triangle_blah)-vert_num-1
                last = 1;
            end

            DOF_mat1(row2:row2+1, col+1) = [triangle_blah(i,7); triangle_blah(i,13)];
            DOF_mat2(row2:row2+1, col+1) = [triangle_blah(i,9-last); triangle_blah(i,15-last)];
        end
        if cyl_def.firstNode ~= "endCap"
            DOF_mat1(:,~any(DOF_mat1,1)) = []; % Remove zero columns
            DOF_mat2(:,~any(DOF_mat2,1)) = [];
        end

    end
    DOF_mat1 = [circshift(DOF_mat1(:,1:end-1), [0 1]), DOF_mat1(:,end)]; % Swap columns, so that DOFs are ascending from column 1
    DOF_mat2 = [circshift(DOF_mat2(:,1:end-1), [0 1]), DOF_mat2(:,end)];
    DOF_mat3 = [circshift(DOF_mat3(:,1:end-1), [0 1]), DOF_mat3(:,end)];
end



temp1        = nonzeros(DOF_mat1); % Create temporary column vector
temp2        = nonzeros(DOF_mat2);
temp3        = nonzeros(DOF_mat3);
% Create vector containing nodes of each edge:
edge_nodes_1 = mesh_data.edges(dof_data.dofs_to_edges(temp1(1:2:end,1)),:); % Edges on contour
edge_nodes_2 = mesh_data.edges(dof_data.dofs_to_edges(temp2(1:2:end,1)),:); % First diagonal
edge_nodes_3 = mesh_data.edges(dof_data.dofs_to_edges(temp3(1:2:end,1)),:); % Second diagonal
clear temp1 temp2 temp3

edge_vecs_1  = mesh_data.node_coords(edge_nodes_1(:,1),:)-mesh_data.node_coords(edge_nodes_1(:,2),:); % Vector of edges
edge_vecs_2  = mesh_data.node_coords(edge_nodes_2(:,1),:)-mesh_data.node_coords(edge_nodes_2(:,2),:);
edge_vecs_3  = mesh_data.node_coords(edge_nodes_3(:,1),:)-mesh_data.node_coords(edge_nodes_3(:,2),:);

 lim1 = length(edge_vecs_1); 
 lim2 = length(edge_vecs_2);
 lim3 = length(edge_vecs_3);
%  if twoEndcaps
% %      lim1 = length(edge_vecs_1) - numVertices;
%  elseif oneEndcap
% %      lim1 = length(edge_vecs_1) - numVertices;
%      if cyl_def.firstNode == "endCap"
% %          lim2 =  length(edge_vecs_2) - numVertices;
%      else
% %          lim3 =  length(edge_vecs_3) - numVertices;
%      end
%  end

theta_1      =  abs(90 - acosd(dot(edge_vecs_1(1:lim2,:),edge_vecs_2(1:lim2,:),2)./(vecnorm(edge_vecs_1(1:lim2,:),2,2).*vecnorm(edge_vecs_2(1:lim2,:),2,2))));
theta_2      =  abs(90 - acosd(dot(edge_vecs_1(1:lim3,:),edge_vecs_3(1:lim3,:),2)./(vecnorm(edge_vecs_1(1:lim3,:),2,2).*vecnorm(edge_vecs_3(1:lim3,:),2,2))));
% quiver3(mesh_data.node_coords(edge_nodes_2(:,1),1),mesh_data.node_coords(edge_nodes_2(:,1),2),mesh_data.node_coords(edge_nodes_2(:,1),3),edge_vecs_2(:,1),edge_vecs_2(:,2),edge_vecs_2(:,3))


for MBF_num = 1:3
    
    B1(1:2,:) = [MBF_mat(edge_nodes_1(:,1),1+MBF_num),MBF_mat(edge_nodes_1(:,2),1+MBF_num)]';
    B2(1:2,:) = [zeros(lim2,1),MBF_mat(edge_nodes_2(:,2),1+MBF_num)]';
    B3(1:2,:) = [MBF_mat(edge_nodes_3(:,1),1+MBF_num),zeros(lim3,1)]';
    B2 = B2(:,:) .* [sind(theta_1)';sind(theta_1)'];
    B3 = B3(:,:) .* [sind(theta_2)';sind(theta_2)'];
    
    for i = 1:lim1
        X1(:,i) = Rho(:,:,i)\B1(:,i);
    end
    for i = 1:lim2 % The diagonal sides have potentially one MBF less if there are endcaps/connections
        X2(:,i) = Rho2(:,:,i)\B2(:,i);
    end
    for i = 1:lim3
        X3(:,i) = Rho2(:,:,i)\B3(:,i);
    end
    
    col_iter = 1;
    node2 = 1;
    node3 = 1;
    for MBF_node = 1:numNodes_new
        if (connection) && MBF_node == numNodes_new
            col_iter = col_iter + 3;
        end
        col_index = col_iter + (MBF_num-1);
        col_iter = col_iter + numMBF;
        % Skip the zero columns
        if DOF_mat1(1,MBF_node) ~= 0
            U_Mat(DOF_mat1(1:2:end,MBF_node),col_index) = X1(1,(numVertices*(MBF_node-1))+1:(numVertices*MBF_node)); % RWG
            U_Mat(DOF_mat1(2:2:end,MBF_node),col_index) = X1(2,(numVertices*(MBF_node-1))+1:(numVertices*MBF_node)); % Linear
        end
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
