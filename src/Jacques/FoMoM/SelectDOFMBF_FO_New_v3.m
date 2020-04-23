function [U_Mat] = SelectDOFMBF_FO_New_v3(mesh_data, dof_data, numVertices ,numMBF, numNodes, triangle_blah, endCap, U_Mat)

% -------------------------------------------------------------------------
% Init
% -------------------------------------------------------------------------
phi = 360/numVertices;
numMBFNodes = (numNodes+2)*numVertices;
DOF_mat = zeros(numVertices*2,numNodes+1);

if endCap
    endCapExclude = (2*numVertices);
else
    endCapExclude = 0;
end


% -------------------------------------------------------------------------
% Set up MBF matrix
% -------------------------------------------------------------------------

sin_mat = sind(phi*(0:(numMBFNodes-1)));
cos_mat = cosd(phi*(0:(numMBFNodes-1)));
contour_nodes = 1:max(max(triangle_blah(:,1:3)));
% contour_nodes = (triangle_blah(1:2:end,3)); % All the nodes associated with the analytical MBF
MBF_mat = [contour_nodes',ones(numMBFNodes,1),sin_mat',cos_mat'];

row = -3;
col = 1;

% Fill each column of matrix with DOFs for each MBF
for i = 1:2:length(triangle_blah)-endCapExclude% Every odd row
    row = row + 4;
    
    DOF_mat(row:row+3,col) = [triangle_blah(i,8);triangle_blah(i,14);triangle_blah(i,7);triangle_blah(i,13)];
    sign1 = triangle_blah(i,5); % Should always be -1
    sign2 = triangle_blah(i,4); % Should always be 1
%     sign3 = 0.5*triangle_blah(i,11);
%     sign4 = 0.5*triangle_blah(i,10);

    Rho(:,:,i) = [1,1;1,-1];
    Rho(:,:,i+1) = [1,1;1,-1];
    if sign1*sign2 == 1
        Rho(:,:,i) = [-1,1;-1,-1];
    end


    if (row == (4*numVertices)-3)
        DOF_mat(row:row+3,col) = [triangle_blah(i,7);triangle_blah(i,13);triangle_blah(i,8);triangle_blah(i,14)];
         Rho(:,:,i+1) = [1,-1;1,1];
        row = -3;
        col = col + 1;
    end
    
end

temp = nonzeros(DOF_mat);
edge_nodes = mesh_data.edges(dof_data.dofs_to_edges(temp(1:2:end,1)),:); % Edges on contour
edge_vecs  = mesh_data.node_coords(edge_nodes(:,1),:)-mesh_data.node_coords(edge_nodes(:,2),:);

% Angles between straight and diagonal edges
theta      =  abs(90 - acosd(dot(edge_vecs(1:2:end,:),edge_vecs(2:2:end,:),2)./(vecnorm(edge_vecs(1:2:end,:),2,2).*vecnorm(edge_vecs(2:2:end,:),2,2))));
% theta      = [theta,theta];

% Retrieve the 1/sin/cos value at the corresponding node
B_const(1:2,:) = [MBF_mat(edge_nodes(:,1),2),MBF_mat(edge_nodes(:,2),2)]';
B_sin(1:2,:)   = [MBF_mat(edge_nodes(:,1),3),MBF_mat(edge_nodes(:,2),3)]';
B_cos(1:2,:)   = [MBF_mat(edge_nodes(:,1),4),MBF_mat(edge_nodes(:,2),4)]';
% Multiply all diagonal edges
B_const(:,2:2:end) = B_const(:,2:2:end) .* [sind(theta)';sind(theta)'];
B_sin(:,2:2:end) = B_sin(:,2:2:end) .* [sind(theta)';sind(theta)'];
B_cos(:,2:2:end) = B_cos(:,2:2:end) .*[sind(theta)';sind(theta)'];


for MBF_num = 1:3
    X = zeros(2,length(edge_nodes));
    switch MBF_num
        case 1
            B = B_const;
        case 2
            B = B_sin;
        case 3
            B = B_cos;
    end

    
    for i = 1:length(edge_nodes) % For every edge
        X(:,i) = Rho(:,:,i)\B(:,i);
    end
    X(2,1:2:end) = 0; % Every second side has a constant value so its linear component should be zero

    col_iter = numMBF-2;
    if col_iter<1
        col_iter = 1;
    end
    for MBF_node = 1:numNodes+1
        col_index = col_iter + (MBF_num-1);
        col_iter = col_iter + numMBF;
        if (numMBF>1) && (ceil(col_index/numMBF)) > numNodes
            col_index = col_index - 3;
        else
            %         col_index = col_index + numMBF;
        end
        
        % x_ind1 = (2*numVertices*(MBF_node-1))+1
        % x_ind2 = (2*numVertices*MBF_node)
        U_Mat(DOF_mat(1:2:end,MBF_node),col_index) = X(1,(2*numVertices*(MBF_node-1))+1:(2*numVertices*MBF_node)); % RWG
        U_Mat(DOF_mat(2:2:end,MBF_node),col_index) = X(2,(2*numVertices*(MBF_node-1))+1:(2*numVertices*MBF_node)); % Linear
    end
    
end