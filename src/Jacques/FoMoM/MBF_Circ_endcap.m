function [U_Mat, DOF_mat, theta1, theta2] = MBF_Circ_endcap(mesh_data, dof_data, numVertices ,numMBF, numNodes, triangle_blah,cyl_def, U_Mat)
% Currently, endcap functionality only works with even number of vertices

% numNodes = numNodes+2;% For coupling
extra = 0;
oneEndcap = (cyl_def.firstNode == "endCap" && cyl_def.lastNode ~= "endCap") || (cyl_def.firstNode ~= "endCap" && cyl_def.lastNode == "endCap"); % Only one is an endcap
twoEndcaps = cyl_def.firstNode == "endCap" && cyl_def.lastNode == "endCap";
connection = cyl_def.firstNode == "conn" && cyl_def.lastNode == "conn";

theta1 = [];
theta2 = [];
% -------------------------------------------------------------------------
% Init
% -------------------------------------------------------------------------
% if cyl_def.firstNode == "conn" || cyl_def.lastNode == "conn"
%    numMBFNodes = (numNodes+3)*numVertices; 
% else
%     numMBFNodes = (numNodes+2)*numVertices;
% end
numMBFNodes = ((numNodes+2)*numVertices) + length(cyl_def.plate_polygon_nodes)+ length(cyl_def.plate_polygon_nodes_end);
phi = 360/numVertices;

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
%         connectionExclude = numVertices; % This means the routine goes an extra 2*numVertices into triangle_blah
        connectionExclude = -numVertices;
%         endCapExclude = numVertices*2;
%         numNodes_new = numNodes + 4;
     end

    
    
elseif twoEndcaps
    extra = 1;
    endCapExclude = (2*numVertices);
    connectionExclude = 0;
else
    connectionExclude = 0;
    if cyl_def.firstNode == "conn" || cyl_def.lastNode == "conn"
        connectionExclude = 0;
    end
    endCapExclude = 0;
    
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

% -------------------------------------------------------------------------
% Set up MBF matrix
% -------------------------------------------------------------------------

sin_mat = circshift(sind(phi*(0:(numMBFNodes-1))),2);
cos_mat = circshift(cosd(phi*(0:(numMBFNodes-1))),2);
 contour_nodes = (1:numMBFNodes)';
ones_mat = (ones(numMBFNodes,1));

if connection
    contour_nodes = (1:(numNodes+2)*numVertices)';
    contour_nodes = [contour_nodes; cyl_def.last_element_val+ cyl_def.plate_polygon_nodes;  cyl_def.last_element_val + cyl_def.plate_polygon_nodes_end];
elseif cyl_def.firstNode == "conn"
    %    contour_nodes(end-numVertices+1:end) = ((numNodes+2)*numVertices) + cyl_def.plate_polygon_nodes ; % Since the polygon nodes are not sequential anymore
     contour_nodes(end-numVertices+1:end) =   ((numNodes+1)*numVertices) + 1 + cyl_def.plate_polygon_nodes ;
%     contour_nodes = [contour_nodes; cyl_def.last_element_val + cyl_def.plate_polygon_nodes];
elseif cyl_def.lastNode == "conn"
    contour_nodes = [contour_nodes; cyl_def.last_element_val + cyl_def.plate_polygon_nodes_end];
end
     
% Contour_nodes can now be used to set the values through linear indexing:
MBF_mat = zeros(numMBFNodes*2, 3); % Don't really care how big this matrix is, as long as the rows correspond to the nodes.
MBF_mat(contour_nodes, :) = [ones_mat,sin_mat',cos_mat'];


% contour_nodes = (triangle_blah(1:2:end,3)); % All the nodes associated with the analytical MBF
% MBF_mat = [contour_nodes',ones_mat,sin_mat',cos_mat'];
if oneEndcap || twoEndcaps
    temp = MBF_mat(:,2) + MBF_mat(:,3);
    [~, maxNode] = max(temp,[],1, 'linear'); % Sin(x) + cos(x)
    if ~twoEndcaps
         MBF_mat((((numNodes+1)*numVertices)+1),:) = [0,0,0];
    else
        MBF_mat((((numNodes+2)*numVertices)+1),:) = [0,0,0]; % The value at the centre vertex, first endcap
        MBF_mat((((numNodes+2)*numVertices)+2),:) = [0,0, 0]; % Second endcap
    end
end


row = -3;
col = 1;
Rho = [1,1;1,-1];
Rho = repmat(Rho, 1,1,length(triangle_blah));

len_tri_mat = length(triangle_blah)-endCapExclude-cyl_def.num_plate_nodes-connectionExclude;

% Fill each column of matrix with DOFs for each MBF
for i = 1:2:len_tri_mat% Every odd row
    row = row + 4;

    sign1 = triangle_blah(i,5); % Should always be -1
    sign2 = triangle_blah(i,4); % Should always be 1
    
    if i <= numVertices*2 && cyl_def.firstNode == "conn"
        DOF_mat(row:row+3,col) = [triangle_blah(i,8);triangle_blah(i,14);triangle_blah(i,9);triangle_blah(i,15)];
    elseif i > len_tri_mat - (2*numVertices) && cyl_def.lastNode == "conn"
        DOF_mat(row:row+3,col) = [triangle_blah(i,8);triangle_blah(i,14);triangle_blah(i,9);triangle_blah(i,15)];
    else
        DOF_mat(row:row+3,col) = [triangle_blah(i,8);triangle_blah(i,14);triangle_blah(i,7);triangle_blah(i,13)];

    end
    
        if sign1*sign2 == 1
            Rho(:,:,i) = [-1,1;-1,-1];
        end
    if (row == (4*numVertices)-3)
        if i <= numVertices*2 && cyl_def.firstNode == "conn"
            DOF_mat(row:row+3,col) = [triangle_blah(i,9);triangle_blah(i,15);triangle_blah(i,8);triangle_blah(i,14)];
        elseif i > len_tri_mat - (2*numVertices) && cyl_def.lastNode == "conn"
            DOF_mat(row:row+3,col) = [triangle_blah(i,9);triangle_blah(i,15);triangle_blah(i,8);triangle_blah(i,14)];
        else
            DOF_mat(row:row+3,col) = [triangle_blah(i,7);triangle_blah(i,13);triangle_blah(i,8);triangle_blah(i,14)];
            Rho(:,:,i+1) = [1,-1;1,1];
        end
        row = -3;
        col = col + 1;
    end
end
if cyl_def.firstNode == "conn" 
    DOF_mat(:,1) = sort(DOF_mat(:,1));
end
if cyl_def.lastNode == "conn" 
    DOF_mat(:,end) = sort(DOF_mat(:,end));
end

ind = 0;
row = -1;
if oneEndcap || twoEndcaps
    if cyl_def.firstNode == "endCap"
        for i = length(triangle_blah)-endCapExclude+1-cyl_def.num_plate_nodes:1:length(triangle_blah)-(numVertices*(cyl_def.lastNode == "endCap")-cyl_def.num_plate_nodes)% Every odd row, only first endcap
            row = row + 2;
            ind = ind + 1;
%            Add the edges on the end cap
            DOF_mat(row:row+1,col) = [triangle_blah(i,7);triangle_blah(i,13)];
            if i == length(triangle_blah)-(numVertices*(cyl_def.lastNode == "endCap")-cyl_def.num_plate_nodes)
                Rho(:,:,i) = [-1,1;-1,-1];
                DOF_mat(row:row+1,col) = [triangle_blah(i,8);triangle_blah(i,14)];
                col = col+1;
            end
        end
    end
    row = -1;
     DOF_mat(:,~any(DOF_mat,1)) = []; % Remove zero columns
    if cyl_def.lastNode == "endCap" %&& cyl_def.firstNode ~= "conn"
        if cyl_def.num_plate_nodes > 0 
            cyl_def.num_plate_nodes = cyl_def.num_plate_nodes -1;
        end
        for i = length(triangle_blah)-cyl_def.num_plate_nodes:-1:length(triangle_blah)-numVertices+1-cyl_def.num_plate_nodes% Second endcap
            row = row + 2; 
            
            if i == length(triangle_blah)-numVertices+1-cyl_def.num_plate_nodes
                Rho(:,:,i) = [-1,-1;-1,1];
                DOF_mat(row:row+1,col) = [triangle_blah(i,7);triangle_blah(i,13)];%;triangle_blah(i+1,8);triangle_blah(i+1,14)]; 
            else
                DOF_mat(row:row+1,col) = [triangle_blah(i,8);triangle_blah(i,14)];  
                if i == length(triangle_blah)-cyl_def.num_plate_nodes
                     Rho(:,:,i) = [1,-1;1,1];
                end
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
    theta(end-round((endCapExclude/2))+1:end,:) = []; % TODO
end
lim = 1:length(edge_nodes);

% Retrieve the 1/sin/cos value at the corresponding node
B_const(1:2,:) = [MBF_mat(edge_nodes(lim,1),1),MBF_mat(edge_nodes(lim,2),1)]';
B_sin(1:2,:)   = [MBF_mat(edge_nodes(lim,1),2),MBF_mat(edge_nodes(lim,2),2)]';
B_cos(1:2,:)   = [MBF_mat(edge_nodes(lim,1),3),MBF_mat(edge_nodes(lim,2),3)]';


% Multiply all diagonal edges
B_const(:,2:2:end-endCapExclude) = B_const(:,2:2:end-endCapExclude) .* [sind(theta)';sind(theta)'];
B_sin(:,2:2:end-endCapExclude) = B_sin(:,2:2:end-endCapExclude) .* [sind(theta)';sind(theta)'];
B_cos(:,2:2:end-endCapExclude) = B_cos(:,2:2:end-endCapExclude) .*[sind(theta)';sind(theta)'];

if oneEndcap || twoEndcaps
    extra = extra + 1;
end


 DOF_mat(:,~any(DOF_mat,1)) = []; % Remove zero columns
        
for MBF_num =1:3
    X = zeros(2,length(edge_nodes));
    switch MBF_num
        case 1
            B = B_const;
        case 2
            B = B_sin;
        case 3
            B = B_cos;
    end

    
    if cyl_def.lastNode == "endCap" % The last endcap's nodes sit on the wrong side for some reason (edge_nodes defined like that)
       B([1 2],end-numVertices+1:end) =   B([2 1],end-numVertices+1:end);
%        if MBF_num == 1 % Fourier mode 0 has phi BC = 0 (at p=0)
%            B(1,end-endCapExclude+1:end) = 0;
%            B(2,end-endCapExclude+1:end) = 0;
%        end
    end
    
    for i = lim % For every edge
        X(:,i) = Rho(:,:,i)\B(:,i);
    end

    col_iter = numMBF-2;
    if col_iter<1
        col_iter = 1;
    end
    x_ind_1 = 0;
%     for MBF_node = 4:4
    for MBF_node = 1:numNodes+1+extra % numNodes is still without extra 2, TODO: add one if this MBF is added to endcap
        col_index = col_iter + (MBF_num-1);
        col_iter = col_iter + numMBF;
        
        fillDofs_1 = nonzeros(DOF_mat(1:2:end, MBF_node));
        fillDofs_2 = nonzeros(DOF_mat(2:2:end, MBF_node));
        
        x_ind_1 = linspace(x_ind_1(end)+1, x_ind_1(end) + length(fillDofs_1),length(fillDofs_1));
        
        U_Mat(fillDofs_1,col_index) = X(1,x_ind_1); % RWG
        U_Mat(fillDofs_2,col_index) = X(2,x_ind_1); % Linear
        
    end
end