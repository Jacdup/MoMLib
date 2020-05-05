function [U_Mat, DOF_mat] = SelectDOFMBF_FO_New_v3(mesh_data, dof_data, numVertices ,numMBF, numNodes, triangle_blah, endCap, U_Mat)
% Currently, endcap functionality only works with even number of vertices


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
contour_nodes = (1:numMBFNodes);
ones_mat = (ones(numMBFNodes,1));

% contour_nodes = 1:max(max(triangle_blah(:,1:3)));
% contour_nodes = (triangle_blah(1:2:end,3)); % All the nodes associated with the analytical MBF
MBF_mat = [contour_nodes',ones_mat,sin_mat',cos_mat'];
if endCap == 1
    [maxVal, maxNode] = max(MBF_mat(:,3) + MBF_mat(:,4),[],1, 'linear'); % Sin(x) + cos(x)
    maxVal = 0;
    MBF_mat(numMBFNodes+1,:) = [(numMBFNodes+1),0,maxVal,maxVal];
    MBF_mat(numMBFNodes+2,:) = [(numMBFNodes+2),0,maxVal,maxVal];
    maxNodes(1,1:2) = [maxNode,numMBFNodes-numVertices+maxNode]; % Nodes where maximum current flows over endcap
    
end




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
next = 0;
if endCap == 1
    
    for i = length(triangle_blah)-endCapExclude+1:2:length(triangle_blah)-numVertices% Every odd row, only first endcap
        row = row + 4;
        
        sign1 = triangle_blah(i,5); % Should always be -1
        sign2 = triangle_blah(i,4); % Should always be 1
        sign3 = triangle_blah(i,10); % Should always be -1
        sign4 = triangle_blah(i,11); % Should always be -1
        
        Rho(:,:,i) = [1,1;1,-1];
        Rho(:,:,i+1) = [1,1;1,-1];
        if sign1*sign2 == 1
            Rho(:,:,i) = [-1,1;-1,-1];
        end
        
        %         if sign3 < 0
        %             Rho(:,:,i) = [1,-1;1,1];
        %         end
        %         if sign4 < 0
        %             Rho(:,:,i+1) = [1,-1;1,1];
        %         end
        
        %         if next == 1
        %             Rho(:,:,i) = [1,-1;1,1];
        %             Rho(:,:,i+1) = [1,-1;1,1];
        %             if sign3 < 0
        %                 Rho(:,:,i) = [1,1;1,-1];
        %             end
        %             if sign4 < 0
        %                 Rho(:,:,i+1) = [1,1;1,-1];
        %             end
        %         end
        
        DOF_mat(row:row+3,col) = [triangle_blah(i,8);triangle_blah(i,14);triangle_blah(i,7);triangle_blah(i,13)];
        
        if (row == (2*numVertices)-3)
            DOF_mat(row:row+3,col) = [triangle_blah(i,7);triangle_blah(i,13);triangle_blah(i,8);triangle_blah(i,14)];
            next = 1;
            %             if sign3 < 0
            %                 Rho(:,:,i+1) = [1,-1;1,1];
            %             end
            %             if sign4 < 0
            %                 Rho(:,:,i) = [1,-1;1,1];
            %             end
            %             Rho(:,:,i+1) = [1,-1;1,1];
            %             if sign3*sign4 < 0
            %                 Rho(:,:,i) = [1,-1;1,1];
            %                 Rho(:,:,i+1) = [1,-1;1,1];
            %             end
            %             if next == 1
            %                 if sign3 < 0
            %                     Rho(:,:,i) = [1,-1;1,1];
            %                 end
            %                 if sign4 < 0
            %                     Rho(:,:,i+1) = [1,-1;1,1];
            %                 end
            %             end
            
            row = -3;
            col = col+1;
        end
    end
    Rho(:,:,i+2) = [1,-1;1,1];
    for i = length(triangle_blah)-numVertices+2:2:length(triangle_blah)% Second endcap
        row = row + 4;
        sign1 = triangle_blah(i,5); % Should always be -1
        sign2 = triangle_blah(i,4); % Should always be 1
        sign3 = triangle_blah(i,10); % Should always be -1
        sign4 = triangle_blah(i,11); % Should always be -1
        DOF_mat(row:row+3,col) = [triangle_blah(i,8);triangle_blah(i,14);triangle_blah(i,7);triangle_blah(i,13)];
        Rho(:,:,i) = [1,-1;1,1];
        Rho(:,:,i+1) = [1,-1;1,1];
        if sign3 < 0
            Rho(:,:,i) = [1,1;1,-1];
        end
        if sign4 < 0
            Rho(:,:,i+1) = [1,1;1,-1];
        end
        if (row == (2*numVertices)-3)
            DOF_mat(row:row+3,col) = [triangle_blah(i,7);triangle_blah(i,13);triangle_blah(i,8);triangle_blah(i,14)];
            next = 1;
            if sign1*sign2 == 1
                Rho(:,:,i+1) = [-1,-1;-1,1];
                if sign3 < 0
                    Rho(:,:,i+1) = [-1,1;-1,-1];
                end
                if sign4 < 0
                    Rho(:,:,i) = [-1,1;-1,-1];
                end
            else
                if sign3 < 0
                    Rho(:,:,i+1) = [1,1;1,-1];
                end
                if sign4 < 0
                    Rho(:,:,i) = [1,1;1,-1];
                end
            end

            row = -3;
            col = col+1;
        end
        
    end
    
end
% Rho(:,:,end-1) = [-1,-1;-1,1];
% Rho(:,:,end)   = [-1,-1;-1,1];




temp = nonzeros(DOF_mat);
edge_nodes = mesh_data.edges(dof_data.dofs_to_edges(temp(1:2:end,1)),:); % Edges on contour
edge_vecs  = mesh_data.node_coords(edge_nodes(:,1),:)-mesh_data.node_coords(edge_nodes(:,2),:);

% edge_nodes(end-(2*numVertices)+1:end,:) = []; % With endcap mesh, but without the circumferential MBF
% edge_vecs(end-(2*numVertices)+1:end,:) = [];

% Angles between straight and diagonal edges
theta      =  abs(90 - acosd(dot(edge_vecs(1:2:end,:),edge_vecs(2:2:end,:),2)./(vecnorm(edge_vecs(1:2:end,:),2,2).*vecnorm(edge_vecs(2:2:end,:),2,2))));
extra = 0;
if endCap == 1
    extra = 1;
    first_endcap = length(edge_nodes)-(2*numVertices)+1:length(edge_nodes)-numVertices;
    second_endcap = length(edge_nodes)-(numVertices)+1:length(edge_nodes);
    
    % find endcap edge connected to maxNode
%     maxEdge1 = edge_nodes(end-(2*numVertices)+1:end-numVertices,:) == maxNodes(1,1); % First endcap
%     maxEdge2 = edge_nodes(end-(numVertices)+1:end,:) == maxNodes(1,2); % Second endcap
    maxEdge1 = edge_nodes(:,:) == maxNodes(1,1); % First endcap
    maxEdge2 = edge_nodes(:,:) == maxNodes(1,2); % Second endcap
    
    max_edge_vec_1 = mesh_data.node_coords(edge_nodes(maxEdge1,1),:)-mesh_data.node_coords(edge_nodes(maxEdge1,2),:);
    max_edge_vec_1 = max_edge_vec_1(3, :); % Only want endcap edge
    max_edge_vec_2 = mesh_data.node_coords(edge_nodes(maxEdge2(:,1),1),:)-mesh_data.node_coords(edge_nodes(maxEdge2(:,1),2),:);
%     max_edge_vec_2 = max_edge_vec_2(3, :); % Only want endcap edge
%     new = repmat(max_edge_vec_1,[numVertices 1]);
    theta1 = abs(90 - acosd(dot(repmat(max_edge_vec_1,[numVertices 1]),edge_vecs(first_endcap,:),2)./(vecnorm(max_edge_vec_1,2,2)*vecnorm(edge_vecs(first_endcap,:),2,2))));
    theta2 = abs(90 -acosd(dot(repmat(max_edge_vec_2,[numVertices 1]),edge_vecs(second_endcap,:),2)./(vecnorm(max_edge_vec_2,2,2)*vecnorm(edge_vecs(second_endcap,:),2,2))));
    
    theta(end-numVertices+1:end,:) = 90; % TODO: dot endcaps with normal component of max
    
    % Angle between every endcap edge and first edge in edge_vecs
end
% theta(numMBFNodes+1: numMBFNodes+1+endCapExclude,:) = 0;
% theta      = [theta,theta];

% Retrieve the 1/sin/cos value at the corresponding node
B_const(1:2,:) = [MBF_mat(edge_nodes(:,1),2),MBF_mat(edge_nodes(:,2),2)]';
B_sin(1:2,:)   = [MBF_mat(edge_nodes(:,1),3),MBF_mat(edge_nodes(:,2),3)]';
B_cos(1:2,:)   = [MBF_mat(edge_nodes(:,1),4),MBF_mat(edge_nodes(:,2),4)]';
% Multiply all diagonal edges  
B_const(:,2:2:end) = B_const(:,2:2:end) .* [sind(theta)';sind(theta)'];
B_sin(:,2:2:end) = B_sin(:,2:2:end) .* [sind(theta)';sind(theta)'];
B_cos(:,2:2:end) = B_cos(:,2:2:end) .*[sind(theta)';sind(theta)'];

if endCap == 1
%    B_const(:,first_endcap) = B_const(:,first_endcap) .* sind(theta1)';
%     B_const(:,second_endcap) = B_const(:,second_endcap) .* sind(theta2)';
   B_const(:,first_endcap) = 0;
    B_const(:,second_endcap) = 0;
       B_sin(:,first_endcap) = B_sin(:,first_endcap) .* sind(theta1)';
    B_sin(:,second_endcap) = B_sin(:,second_endcap) .* sind(theta2)';
       B_cos(:,first_endcap) = B_cos(:,first_endcap) .* sind(theta1)';
    B_cos(:,second_endcap) = B_cos(:,second_endcap) .* sind(theta2)';
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

    
    for i = 1:length(edge_nodes) % For every edge
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
            col_index = col_index - 3;
        else
            %         col_index = col_index + numMBF;
        end
        
        % x_ind1 = (2*numVertices*(MBF_node-1))+1
        % x_ind2 = (2*numVertices*MBF_node)
        if MBF_node > numNodes + 1 % Last two
            col_index = col_index + 3;
            U_Mat(DOF_mat(1:2:end/2,MBF_node),col_index) = X(1,(2*numVertices*(MBF_node-1))+1:(2*numVertices*MBF_node)-numVertices); % RWG
            U_Mat(DOF_mat(2:2:end/2,MBF_node),col_index) = X(2,(2*numVertices*(MBF_node-1))+1:(2*numVertices*MBF_node)-numVertices); % Linear
            col_index = col_index + 3;
            U_Mat(DOF_mat(1:2:end/2,MBF_node+1),col_index) = X(1,(2*numVertices*(MBF_node-1))+1+numVertices:(2*numVertices*MBF_node)); % RWG
            U_Mat(DOF_mat(2:2:end/2,MBF_node+1),col_index) = X(2,(2*numVertices*(MBF_node-1))+1+numVertices:(2*numVertices*MBF_node)); % Linear
        else
            U_Mat(DOF_mat(1:2:end,MBF_node),col_index) = X(1,(2*numVertices*(MBF_node-1))+1:(2*numVertices*MBF_node)); % RWG
            U_Mat(DOF_mat(2:2:end,MBF_node),col_index) = X(2,(2*numVertices*(MBF_node-1))+1:(2*numVertices*MBF_node)); % Linear
        end
        
        

    end
    
end