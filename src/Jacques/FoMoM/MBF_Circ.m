function [U_Mat, DOF_mat, theta1, theta2] = MBF_Circ(mesh_data, dof_data, numVertices ,numMBF, numNodes, triangle_blah, endCap, connection,cyl_def, U_Mat)
% Currently, endcap functionality only works with even number of vertices

% numNodes = numNodes+2;% For coupling
extra = 0;
oneEndcap = (cyl_def.firstNode == "endCap" && cyl_def.lastNode ~= "endCap") || (cyl_def.firstNode ~= "endCap" && cyl_def.lastNode == "endCap"); % Only one is an endcap
twoEndcaps = cyl_def.firstNode == "endCap" && cyl_def.lastNode == "endCap";

theta1 = [];
theta2 = [];
% -------------------------------------------------------------------------
% Init
% -------------------------------------------------------------------------
phi = 360/numVertices;
numMBFNodes = (numNodes+2)*numVertices;
DOF_mat = zeros(numVertices*2,numNodes+1);
% theta1 = 
% if connection 
%     endCapExclude = 4*numVertices;
% else
if oneEndcap % TODO: connection functionality
    endCapExclude = (numVertices);
    connectionExclude = 0;
    
     if cyl_def.firstNode == "conn" || cyl_def.lastNode == "conn"
%         connCapExclude = vert_num;
        connectionExclude = numVertices;
%         endCapExclude = numVertices*2;
%         numNodes_new = numNodes + 4;
     end

    
    
elseif twoEndcaps
    endCapExclude = (2*numVertices);
    connectionExclude = 0;
else
    endCapExclude = 0;
    connectionExclude = 0;
end

% if endCap && (connection==0)
%     endCapExclude = (2*numVertices);
%     connectionExclude = 0;
% elseif connection
%     endCapExclude = 0;
%     connectionExclude = (2*numVertices);
% else
%     endCapExclude = 0;
%     connectionExclude = 0;
% end

EPS = 0.0001;
% -------------------------------------------------------------------------
% Set up MBF matrix
% -------------------------------------------------------------------------

sin_mat = sind(phi*(0:(numMBFNodes-1)));
cos_mat = cosd(phi*(0:(numMBFNodes-1)));
contour_nodes = (1:numMBFNodes);
ones_mat = (ones(numMBFNodes,1));

% contour_nodes = 1:max(max(triangle_blah(:,1:3)));
% contour_nodes = (triangle_blah(1:2:end,3)); % All the nodes associated with the analytical MBF
MBF_mat = [contour_nodes',ones_mat,sin_mat',cos_mat'];
if oneEndcap || twoEndcaps
    temp = MBF_mat(:,3) + MBF_mat(:,4);
    [maxVal, maxNode] = max(temp,[],1, 'linear'); % Sin(x) + cos(x)
    maxVal = 0; % The value at the centre of the endcap
    MBF_mat(numMBFNodes+1,:) = [(numMBFNodes+1),0,maxVal,maxVal];
    MBF_mat(numMBFNodes+2,:) = [(numMBFNodes+2),0,maxVal,maxVal];
    maxNodes(1,1:2) = [maxNode,numMBFNodes-numVertices+maxNode]; % Nodes where maximum current flows over endcap
end


row = -3;
col = 1;
Rho = [1,1;1,-1];
Rho = repmat(Rho, 1,1,length(triangle_blah));

% Fill each column of matrix with DOFs for each MBF
for i = 1:2:length(triangle_blah)-endCapExclude-connectionExclude-cyl_def.num_plate_nodes% Every odd row
    row = row + 4;
    DOF_mat(row:row+3,col) = [triangle_blah(i,8);triangle_blah(i,14);triangle_blah(i,7);triangle_blah(i,13)];
    sign1 = triangle_blah(i,5); % Should always be -1
    sign2 = triangle_blah(i,4); % Should always be 1
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

if oneEndcap || twoEndcaps
    if cyl_def.firstNode == "endCap"
        for i = length(triangle_blah)-endCapExclude+1-cyl_def.num_plate_nodes:2:length(triangle_blah)-(numVertices*(cyl_def.lastNode == "endCap")-cyl_def.num_plate_nodes)% Every odd row, only first endcap
            row = row + 4;
            sign1 = triangle_blah(i,5); % Should always be -1
            sign2 = triangle_blah(i,4); % Should always be 1
            if sign1*sign2 == 1
                Rho(:,:,i) = [-1,1;-1,-1];
            end
            DOF_mat(row:row+3,col) = [triangle_blah(i,8);triangle_blah(i,14);triangle_blah(i,7);triangle_blah(i,13)];
            if (row == (2*numVertices)-3)
                DOF_mat(row:row+3,col) = [triangle_blah(i,7);triangle_blah(i,13);triangle_blah(i,8);triangle_blah(i,14)];
                row = -3;
                col = col+1;
            end
        end
    end
%     Rho(:,:,i+2) = [1,-1;1,1];
    if cyl_def.lastNode == "endCap"
        for i = length(triangle_blah)-numVertices+2-cyl_def.num_plate_nodes:2:length(triangle_blah)-cyl_def.num_plate_nodes% Second endcap
            row = row + 4; 
            DOF_mat(row:row+3,col) = [triangle_blah(i,8);triangle_blah(i,14);triangle_blah(i,7);triangle_blah(i,13)];    
            [Rho(:,:,i), Rho(:,:,i+1)] = getSigns(triangle_blah(i,:));   
            if (row == (2*numVertices)-3)
                DOF_mat(row:row+3,col) = [triangle_blah(i,7);triangle_blah(i,13);triangle_blah(i,8);triangle_blah(i,14)];
                [Rho(:,:,i), Rho(:,:,i+1)] = getSigns(triangle_blah(i,:));          
                row = -3;
                col = col+1;
            end
        end
    end
end


temp = nonzeros(DOF_mat);
edge_nodes = mesh_data.edges(dof_data.dofs_to_edges(temp(1:2:end,1)),:); % Edges on contour
edge_vecs  = mesh_data.node_coords(edge_nodes(:,1),:)-mesh_data.node_coords(edge_nodes(:,2),:);

% edge_nodes(end-(2*numVertices)+1:end,:) = []; % With endcap mesh, but without the circumferential MBF
% edge_vecs(end-(2*numVertices)+1:end,:) = [];

% Angles between straight and diagonal edges
theta      =  abs(90 - acosd(dot(edge_vecs(1:2:end,:),edge_vecs(2:2:end,:),2)./(vecnorm(edge_vecs(1:2:end,:),2,2).*vecnorm(edge_vecs(2:2:end,:),2,2))));

if oneEndcap || twoEndcaps
    extra = 1;
    
    if cyl_def.firstNode == "endCap"
        first_endcap = length(edge_nodes)-(2*numVertices)+1:length(edge_nodes)-numVertices;
%         maxEdge1 = edge_nodes(:,:) == maxNodes(1,1); % First endcap
        maxFirstEndCapNode1 = mesh_data.node_coords(edge_nodes(maxNode,1),:); % Smallest numbers
        maxFirstEndCapNode1 = maxFirstEndCapNode1(1,:); % First node
        maxFirstEndCapNode2 = mesh_data.node_coords(end-1,:);
        max_edge_vec_1 = maxFirstEndCapNode1 - maxFirstEndCapNode2;
         % Angle between every endcap edge and first edge in edge_vecs
        theta1 = abs(90 - acosd(dot(repmat(max_edge_vec_1,[numVertices 1]),edge_vecs(first_endcap,:),2)./(vecnorm(max_edge_vec_1,2,2)*vecnorm(edge_vecs(first_endcap,:),2,2))));
    end
    
    if cyl_def.lastNode == "endCap"
        second_endcap = length(edge_nodes)-(numVertices)+1:length(edge_nodes);
%         maxEdge2 = edge_nodes(:,:) == maxNodes(1,2); % Second endcap
        maxSecondEndCapNode1 = mesh_data.node_coords(edge_nodes(maxNode(:,1),1),:);
%         maxSecondEndCapNode1 = mesh_data.node_coords(edge_nodes(maxEdge2(:,1),1),:);
        maxSecondEndCapNode1 = maxSecondEndCapNode1(end,:); % Last node
        maxSecondEndCapNode2 = mesh_data.node_coords(end,:);
        max_edge_vec_2 = maxSecondEndCapNode1 - maxSecondEndCapNode2;
        %      max_edge_vec_2 = max_edge_vec_1;
         % Angle between every endcap edge and first edge in edge_vecs
        theta2 = abs(90 - acosd(dot(repmat(max_edge_vec_2,[numVertices 1]),edge_vecs(second_endcap,:),2)./(vecnorm(max_edge_vec_2,2,2)*vecnorm(edge_vecs(second_endcap,:),2,2))));
    end
    theta(end-round((endCapExclude/2))+1:end,:) = []; % TODO
end
lim = 1:length(edge_nodes);
if cyl_def.firstNode == "conn"
   lim = 1:length(edge_nodes)- numVertices;
   theta(end-round((endCapExclude/2))+1:end,:) = []; % TODO
end
% Retrieve the 1/sin/cos value at the corresponding node
B_const(1:2,:) = [MBF_mat(edge_nodes(lim,1),2),MBF_mat(edge_nodes(lim,2),2)]';
B_sin(1:2,:)   = [MBF_mat(edge_nodes(lim,1),3),MBF_mat(edge_nodes(lim,2),3)]';
B_cos(1:2,:)   = [MBF_mat(edge_nodes(lim,1),4),MBF_mat(edge_nodes(lim,2),4)]';

% Multiply all diagonal edges
B_const(:,2:2:end-endCapExclude) = B_const(:,2:2:end-endCapExclude) .* [sind(theta)';sind(theta)'];
B_sin(:,2:2:end-endCapExclude) = B_sin(:,2:2:end-endCapExclude) .* [sind(theta)';sind(theta)'];
B_cos(:,2:2:end-endCapExclude) = B_cos(:,2:2:end-endCapExclude) .*[sind(theta)';sind(theta)'];
if oneEndcap || twoEndcaps

    %    B_const(:,first_endcap) = B_const(:,first_endcap) .* sind(theta1)';
    %     B_const(:,second_endcap) = B_const(:,second_endcap) .* sind(theta2)';
    if cyl_def.firstNode == "endCap"
        B_const(:,first_endcap) = 0;
    end
    if cyl_def.lastNode == "endCap"
        B_const(:,second_endcap) = 0;
    end
    
    %     B_sin(:,first_endcap)    = (B_sin(:,first_endcap)  .* (sind(theta1))')  ;
    %     B_sin(:,second_endcap)   = (B_sin(:,second_endcap) .* (sind(theta2))') ;
    %     B_cos(:,first_endcap)    = (B_cos(:,first_endcap)  .* (sind(theta1))') ;
    %     B_cos(:,second_endcap)   = (B_cos(:,second_endcap) .* (sind(theta2))');
end


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

    
    for i = lim % For every edge
        X(:,i) = Rho(:,:,i)\B(:,i);
    end
%     X(2,1:2:end) = 0; % Every second side has a constant value so its linear component should be zero

    col_iter = numMBF-2;
    if col_iter<1
        col_iter = 1;
    end
    for MBF_node = 1:numNodes+1+extra % numNodes is still without extra 2, TODO: add one if this MBF is added to endcap
        col_index = col_iter + (MBF_num-1);
        col_iter = col_iter + numMBF;
        if (numMBF>1) && (ceil(col_index/numMBF)) > numNodes+extra
            if connection
                col_index = col_index + 3;
            else
                col_index = col_index - 3;
            end
        else
%                     col_index = col_index + numMBF;
        end
        
        % x_ind1 = (2*numVertices*(MBF_node-1))+1
        % x_ind2 = (2*numVertices*MBF_node)
        if (oneEndcap || twoEndcaps) && (MBF_node > numNodes + 1) % Last two
            col_index = col_index + 3;
            if cyl_def.firstNode == "endCap"
                U_Mat(DOF_mat(1:2:end/2,MBF_node),col_index) = X(1,(2*numVertices*(MBF_node-1))+1:(2*numVertices*MBF_node)-numVertices); % RWG
                U_Mat(DOF_mat(2:2:end/2,MBF_node),col_index) = X(2,(2*numVertices*(MBF_node-1))+1:(2*numVertices*MBF_node)-numVertices); % Linear
                col_index = col_index + 3;
            end
            if (cyl_def.lastNode == "endCap" && cyl_def.firstNode ~= "endCap")
                U_Mat(DOF_mat(1:2:end/2,MBF_node),col_index) = X(1,(2*numVertices*(MBF_node-1))+1:(2*numVertices*MBF_node)-numVertices); % RWG
                U_Mat(DOF_mat(2:2:end/2,MBF_node),col_index) = X(2,(2*numVertices*(MBF_node-1))+1:(2*numVertices*MBF_node)-numVertices); % Linear
            end
            if twoEndcaps
                U_Mat(DOF_mat(1:2:end/2,MBF_node+1),col_index) = X(1,(2*numVertices*(MBF_node-1))+1+numVertices:(2*numVertices*MBF_node)); % RWG
                U_Mat(DOF_mat(2:2:end/2,MBF_node+1),col_index) = X(2,(2*numVertices*(MBF_node-1))+1+numVertices:(2*numVertices*MBF_node)); % Linear
            end
%                 end
        else
            U_Mat(DOF_mat(1:2:end,MBF_node),col_index) = X(1,(2*numVertices*(MBF_node-1))+1:(2*numVertices*MBF_node)); % RWG
            U_Mat(DOF_mat(2:2:end,MBF_node),col_index) = X(2,(2*numVertices*(MBF_node-1))+1:(2*numVertices*MBF_node)); % Linear
        end
        
    end
end

% 
% if oneEndcap && twoEndcaps
%     U_Mat(:,end-2) = []; % Remove constant columns at endcaps
%     U_Mat(:,end-4) = [];
% end