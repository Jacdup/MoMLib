function [U_Mat, DOF_mat1, DOF_mat2, DOF_mat3] = MBF_Axial_endcap(mesh_data, dof_data, numVertices ,numMBF, numNodes, triangle_blah, cyl_def)
% TODO: make DOF_Mat1/2/3 a single cell matrix, so that it is much easier
% to figure out which columns corresponds/which DOFs are part of the same
% MBF.
% The MBF is agnostic to signs on the basis functions, and only cares about
% if all the DOFS where the MBF is defined have the same orientation

oneEndcap = (cyl_def.firstNode == "endCap" && cyl_def.lastNode ~= "endCap") || (cyl_def.firstNode ~= "endCap" && cyl_def.lastNode == "endCap"); % Only one is an endcap
twoEndcaps = cyl_def.firstNode == "endCap" && cyl_def.lastNode == "endCap";
connection = cyl_def.firstNode == "conn" && cyl_def.lastNode == "conn";
% oneEndcap = 0;
% twoEndcaps = 0;
% -------------------------------------------------------------------------
% Init
% -------------------------------------------------------------------------
phi = 360/numVertices;
vert_num = (2*numVertices)-1;
extra_dof_col = 0;
if twoEndcaps
    % if endCap || connection
    % In the connection case, the setup is exactly the same as the endcap,
    % since it is a special case of the endcap
    numNodes_new = numNodes + 2;
    endCapExclude = (2*numVertices); % exclude last (2*numVertices) from i/row assignment
    connCapExclude = 0;
elseif oneEndcap
    numNodes_new = numNodes + 1;

    if cyl_def.firstNode == "conn" || cyl_def.lastNode == "conn"
%         connCapExclude = vert_num;
        connCapExclude = -vert_num -1;
%         endCapExclude = numVertices;
        endCapExclude = (3*numVertices);
%         numNodes_new = numNodes + 4;
    numNodes_new = numNodes + 2;
%      numNodes_new = numNodes + 1;
    else
        endCapExclude = numVertices;
        connCapExclude = 0;
    end
    
elseif connection
     connCapExclude = -(vert_num*2)-2;
     endCapExclude = (4*numVertices);
     numNodes_new = numNodes + 3;
else
    numNodes_new = numNodes;
    endCapExclude = 0;
    connCapExclude = 0;
end
% temp = 1;
numDofs = size(dof_data.basis_supports,1);
numMBFNodes = numVertices*(numNodes_new) ;
DOF_mat1 = zeros(numVertices*2,numNodes_new);
DOF_mat2 = zeros(numVertices*2,numNodes_new);
DOF_mat3 = zeros(numVertices*2,numNodes_new);
X1       = zeros(2,numMBFNodes);
X2       = zeros(2,numMBFNodes);
X3       = zeros(2,numMBFNodes);
U_Mat = zeros(numDofs, numNodes_new*numMBF);

len_tri_mat = (length(triangle_blah)-vert_num-1-endCapExclude-cyl_def.num_plate_nodes - connCapExclude);

Rho = [1,1;1,-1]; % This is the matrix of RWG and linear components at the two edge nodes
Rho = repmat(Rho,1,1,numMBFNodes);
Rho2 = Rho;
Rho3 = Rho;
% -------------------------------------------------------------------------
% Set up MBF matrix
% -------------------------------------------------------------------------
% MBF matrix is the analytical MBF
% MBF_mat has the value of the MBF at each contour node point
% num_nodes x [nodes,constant,sin,cos]
% numMBFNodes_new = numVertices*(numNodes+4); % For coupling
% numMBFNodes_new = numVertices*(numNodes_new+2) ; % This is only for filling the MBF_mat
numMBFNodes_new = ((numNodes+2)*numVertices) + length(cyl_def.plate_polygon_nodes) + length(cyl_def.plate_polygon_nodes_end);
sin_mat = sind(phi*(0:(numMBFNodes_new-1)));
cos_mat = cosd(phi*(0:(numMBFNodes_new-1)));
ones_mat = (ones(numMBFNodes_new,1));

contour_nodes = (1:numMBFNodes_new)';

if connection
    contour_nodes = (1:(numNodes+2)*numVertices)';
    %     contour_nodes(end-(numVertices*2)+1:end-numVertices) =   ((numNodes+1)*numVertices) + 1 + cyl_def.plate_polygon_nodes ;
    %     contour_nodes(end-(numVertices)+1:end)               =   ((numNodes+1)*numVertices) + 1 + cyl_def.plate_polygon_nodes_end ;
    %     * (cyl_def.lastNode ~= "conn"))
    contour_nodes = [contour_nodes; cyl_def.last_element_val+  cyl_def.plate_polygon_nodes;  cyl_def.last_element_val + cyl_def.plate_polygon_nodes_end];
    
elseif cyl_def.firstNode == "conn"
    %     contour_nodes = [contour_nodes; cyl_def.last_element_val  + cyl_def.plate_polygon_nodes];
    contour_nodes(end-numVertices+1:end) =   ((numNodes+1)*numVertices) + 1 + cyl_def.plate_polygon_nodes ; % Since the polygon nodes are not sequential anymore
elseif cyl_def.lastNode == "conn"
    contour_nodes = [contour_nodes; cyl_def.last_element_val + cyl_def.plate_polygon_nodes_end];
    %     contour_nodes(end-numVertices+1:end) =   ((numNodes+1)*numVertices) + 1 + cyl_def.plate_polygon_nodes_end ;
end
% Contour_nodes can now be used to set the values through linear indexing:
MBF_mat = zeros(numMBFNodes*2, 3); % Don't really care how big this matrix is, as long as the rows correspond to the nodes.
MBF_mat(contour_nodes, :) = [ones_mat,sin_mat',cos_mat'];
% ------------------------------------------------------------------------
% New 05/11/2020
% ------------------------------------------------------------------------
if oneEndcap || twoEndcaps
%     temp = MBF_mat(:,2) + MBF_mat(:,3);
%     [~, maxNode] = max(temp,[],1, 'linear'); % Sin(x) + cos(x)
    if ~twoEndcaps
         MBF_mat((((numNodes+1)*numVertices)+1),:) = [0,1,1]; % End cap vertex always has value of unity (max)
%          MBF_mat((((numNodes+1)*numVertices)+2),:) = [0,1,1];
    else
        MBF_mat((((numNodes+2)*numVertices)+1),:) = [0,1,1]; % The value at the centre vertex, first endcap
        MBF_mat((((numNodes+2)*numVertices)+2),:) = [0,1,1];
    end
%     maxNodes(1,1:2) = [maxNode,numMBFNodes-numVertices+maxNode]; % Nodes where maximum current flows over endcap
end

% MBF_mat = MBF_mat(1:end-36,:); % Temporary for coupling

row = -1;
iter = 1;
col = 1;
linear_row = 0;
rho_ind = 0;

%  [DOF_mat1,DOF_mat2,DOF_mat3] =
%  getAxialDOFs(triangle_blah,numVertices,numNodes, len_tri_mat); % TODO,
%  input stuff, get out DOFs of axial MBF, ordered.

% Select DOFs associated with the MBF
for i = 1:2:len_tri_mat
    linear_row = linear_row + 1;
    rho_ind = rho_ind + 1;
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
    
    if i <= numVertices*2 && cyl_def.firstNode == "conn" % Assign first MBF with new DOFs from generated connection
        % This is NOT the transition from cyl -> plate
        DOF_mat1(row:row+1,col) = [triangle_blah(i+1,9);triangle_blah(i+1,15)];
        DOF_mat2(row:row+1,col) = [triangle_blah(i+1-last,8);triangle_blah(i+1-last,14)];
    end
    if cyl_def.lastNode == "conn" && i >= len_tri_mat - vert_num +1 % Assign last MBF with new DOFs from generated connection
        % This is NOT the transition from cyl -> plate
        DOF_mat3(row:row+1,col) = [triangle_blah(i+vert_num+2-last,8);triangle_blah(i+vert_num+2-last,14)];
    end

    if last == 1
        row = -1;
        linear_row = 0;
        col = col + 1;
    end
end
for i = numVertices:numVertices:length(Rho) % DOFs at the end of each node have switched signs
    Rho2(:,:,i) = Rho2(:,:,i) .* [1,-1;1,-1];
     Rho3(:,:,i) = Rho3(:,:,i) .* [1,-1;1,-1];
    %     Rho(:,:,i) = Rho(:,:,i) .* [1,-1;1,-1];
end

% Temporary for coupling
% DOF_mat1(:,32) = [];
% DOF_mat2(:,32) = [];
% DOF_mat3(:,32) = [];
if connection
    
%     Rho2(:,:,end-numVertices+1:end) = Rho2(:,:,end-numVertices+1:end) .* [-1,-1;-1,-1];
%     Rho3(:,:,end-numVertices+1:end) = Rho3(:,:,end-numVertices+1:end) .* [-1,-1;-1,-1];
    
%     Rho2(:,:,end-numVertices+1:end) = Rho2(:,:,end-numVertices+1:end) .* [1,-1;1,-1];
%     Rho3(:,:,1:numVertices) = Rho3(:,:,1:numVertices) .* [-1,1;-1,1];
end
% -------------------------------------------------------------------------
% -----------------------------Endcap stuff--------------------------------
% -------------------------------------------------------------------------
if oneEndcap || twoEndcaps || connection
    row1 = -1; % Fill up a new column, from row 1
    row2 = -1;
    linear_row = 0;
    if cyl_def.firstNode == "endCap" || cyl_def.firstNode == "conn"
        for i = 1:2:(vert_num) % First endcap, or connection triangles
            linear_row = linear_row + 1;
            row1 = row1+ 2;
            last = 0;
            if (i == vert_num)
                last = 1;
            end

            if cyl_def.firstNode == "conn"
                DOF_mat1(row1:row1+1, col) = [triangle_blah(i,7); triangle_blah(i,13)];
                DOF_mat3(row1:row1+1, col) = [triangle_blah(i+1,8-last); triangle_blah(i+1,14-last)];
                
%                  Rho2(:,:, i) = Rho2(:,:,i) .* [-1,-1;-1,-1];
            if connection
                Rho_index = sub2ind(size(DOF_mat1(1:2:end,:)), linear_row, col);
                Rho3(:,:, Rho_index) = Rho3(:,:,Rho_index) .* [-1,-1;-1,-1]; % This fixes the 'connection' case (2 connections)
            end
%                  
            else
%                  Rho_index = sub2ind(size(DOF_mat1(1:2:end,:)), linear_row, col);
%                  Rho(:,:, Rho_index) = Rho(:,:,Rho_index) .* [-1,-1;-1,-1];
                DOF_mat1(row1:row1+1, col) = [triangle_blah(i,9); triangle_blah(i,15)];
                DOF_mat3(row1:row1+1, col) = [triangle_blah(i,7+last); triangle_blah(i,13+last)];
            end
            
        end
%         DOF_mat1(:,col) = sort(DOF_mat1(:,col));
    end
    
    linear_row = 0;
    if cyl_def.lastNode == "endCap" || cyl_def.lastNode == "conn" 
        for i = (length(triangle_blah)-vert_num-endCapExclude- connCapExclude - cyl_def.num_plate_nodes+1):2:length(triangle_blah)-endCapExclude-connCapExclude-cyl_def.num_plate_nodes % Second endcap, or connection triangles
            row2 = row2 + 2;
            linear_row = linear_row + 1;
            last = 0;
            
            if i ==  length(triangle_blah)-endCapExclude-connCapExclude-cyl_def.num_plate_nodes
                last = 1;
            end
%             Rho_index = sub2ind(size(DOF_mat1(1:2:end,:)), linear_row, col+1);
%             DOF_mat1(row2:row2+1, col+1+extra_dof_col) = [triangle_blah(i,7); triangle_blah(i,13)];
%             if cyl_def.lastNode == "conn" 
            if cyl_def.firstNode == "conn" && cyl_def.lastNode ~= "conn" % Dofs are assigned differently when first node is a connection
               DOF_mat2(row2:row2+1, col+1) = [triangle_blah(i+last,7+last); triangle_blah(i+last,13+last)];
                DOF_mat1(row2:row2+1, col+1) = [triangle_blah(i+1,7); triangle_blah(i+1,13)];
%                 Rho_index = sub2ind(size(DOF_mat1(1:2:end,:)), linear_row, col+1);
                Rho3(2,2,linear_row) = -0.5*triangle_blah(i+last,13+last);
                Rho3(1,2,linear_row) = 0.5*triangle_blah(i+last,13+last);
                Rho3(1,1,linear_row) = 1;
                 Rho3(2,1,linear_row) = 1;
                 Rho2(:,:,linear_row) = Rho3(:,:,linear_row);
%                 Rho3(:,:, linear_row) = Rho3(:,:,linear_row) .* [1,-1;1,-1]; % This fixes the 'connection' case (2 connections)
            elseif cyl_def.lastNode == "conn"
                DOF_mat1(row2:row2+1, col+1) = [triangle_blah(i,7); triangle_blah(i,13)];
                DOF_mat2(row2:row2+1, col+1) = [triangle_blah(i+1-last,8); triangle_blah(i+1-last,14)];
                if connection
                    Rho_index = sub2ind(size(DOF_mat1(1:2:end,:)), linear_row, col+1);
                    Rho(:,:, Rho_index) = Rho(:,:,Rho_index) .* [-1,-1;-1,-1]; % This fixes the 'connection' case (2 connections)
%                     Rho2(:,:, Rho_index) = Rho2(:,:,Rho_index) .* [-1,-1;-1,-1];
%                     Rho3(:,:, Rho_index) = Rho3(:,:,Rho_index) .* [-1,-1;-1,-1];% uncomment this for pretty current
                end
%                 Rho3(:,:, Rho_index) = Rho3(:,:,Rho_index) .* [-1,-1;-1,-1];
%                 Rho2(:,:, Rho_index) = Rho2(:,:,Rho_index) .* [1,-1;1,-1];
            else
                DOF_mat1(row2:row2+1, col+1) = [triangle_blah(i,7); triangle_blah(i,13)];
                DOF_mat2(row2:row2+1, col+1) = [triangle_blah(i,9-last); triangle_blah(i,15-last)];
            end

%              Rho3(:,:, Rho_index) = Rho3(:,:,Rho_index) .* [-1,-1;-1,-1];
%             [Rho(:,:,Rho_index), Rho2(:,:,Rho_index)] = getSigns_new(triangle_blah,7,9-last,i,i);
            
        end
        if cyl_def.firstNode ~= "endCap" && cyl_def.firstNode ~= "conn" % If the first node is hollow
            DOF_mat1(:,~any(DOF_mat1,1)) = []; % Remove zero columns
            DOF_mat2(:,~any(DOF_mat2,1)) = [];
        end
        
    end
end
% -------------------------------------------------------------------------
% This comes straight from the Circ Endcap, to select the edges on the
% endcap
ind = 0;
if oneEndcap || twoEndcaps
    row = -1;
    if cyl_def.firstNode == "endCap"
        if twoEndcaps
            col_ind = 1;
            DOF_mat2 = [circshift(DOF_mat2(:,1:end-1), [0 1]), DOF_mat2(:,end)];
        else
            col_ind = col;
        end
        for i = length(triangle_blah)-endCapExclude+1-cyl_def.num_plate_nodes:1:length(triangle_blah)-(numVertices*(cyl_def.lastNode == "endCap")-cyl_def.num_plate_nodes)% Every odd row, only first endcap
            row = row + 2;
            ind = ind + 1;
%             Rho2(:,:,ind) = [-1,1;-1,-1];
%            Add the edges on the end cap
            DOF_mat2(row:row+1,col_ind) = [triangle_blah(i,7);triangle_blah(i,13)];
            if i == length(triangle_blah)-(numVertices*(cyl_def.lastNode == "endCap")-cyl_def.num_plate_nodes)
                Rho2(:,:,ind) = [-1,1;-1,-1];
                 DOF_mat2(row:row+1,col_ind) = [triangle_blah(i,8);triangle_blah(i,14)];
            end
        end
    end
    col = col+1;
%     DOF_mat2(1:numVertices*2,1) = circshift(DOF_mat2(1:numVertices*2,1),2);
%      DOF_mat2(:,~any(DOF_mat,1)) = []; % Remove zero columns
%     Rho(:,:,i+2) = [1,-1;1,1];
    row = -1;
    ind = length(Rho3) + 1;
    if cyl_def.lastNode == "endCap" %&& cyl_def.firstNode ~= "conn"
        if cyl_def.num_plate_nodes > 0 
            cyl_def.num_plate_nodes = cyl_def.num_plate_nodes -1;
        end
        for i = length(triangle_blah)-cyl_def.num_plate_nodes:-1:length(triangle_blah)-numVertices+1-cyl_def.num_plate_nodes% Second endcap
            row = row + 2; 
            ind = ind - 1;
            if i == length(triangle_blah)-numVertices+1-cyl_def.num_plate_nodes
                Rho3(:,:,ind) = [-1,-1;-1,1];
                DOF_mat3(row:row+1,col) = [triangle_blah(i,7);triangle_blah(i,13)];%;triangle_blah(i+1,8);triangle_blah(i+1,14)]; 
            else
                DOF_mat3(row:row+1,col) = [triangle_blah(i,8);triangle_blah(i,14)];  
                if i == length(triangle_blah)-cyl_def.num_plate_nodes
                     Rho3(:,:,ind) = [1,-1;1,1];
                end
            end

        end
    end
end


% -------------------------------------------------------------------------
if twoEndcaps || (cyl_def.firstNode == "conn" && cyl_def.lastNode ~= "conn")
    DOF_mat1 = [circshift(DOF_mat1(:,1:end-1), [0 1]), DOF_mat1(:,end)]; % Swap columns, so that DOFs are ascending from column 1
    if cyl_def.firstNode == "conn"
         DOF_mat2 = [circshift(DOF_mat2(:,1:end-1), [0 1]), DOF_mat2(:,end)];
    end
    DOF_mat3 = [circshift(DOF_mat3(:,1:end-1), [0 1]), DOF_mat3(:,end)];
    Rho(:,:,1:numVertices) = Rho(:,:,1:numVertices) .* [-1,-1;-1,-1]; % I REALLY don't know why this is suddenly necessary (since 25/05/2020)
elseif cyl_def.firstNode == "endCap"
%       Rho(:,:,1:numVertices) = Rho(:,:,1:numVertices) .* [-1,-1;-1,-1]; 
    Rho(:,:,end-numVertices+1:end) = Rho(:,:,end-numVertices+1:end) .* [-1,-1;-1,-1]; % I REALLY don't know why this is suddenly necessary (since 25/05/2020)
%     Rho2(:,:,1:numVertices) = Rho2(:,:,1:numVertices) .* [-1,-1;-1,-1];
end
% Create cell array that has all the DOFs associated with all the MBFs 
% (One MBF per cell)

if twoEndcaps
%     DOF_mat2( :, ~any(DOF_mat2,1) ) = [];  %remove zero columns
%     DOF_mat3( :, ~any(DOF_mat3,1) ) = [];  %remove zero column
end

for k = 1:size(DOF_mat1,2)
%     DOF_mat{k} = [DOF_mat2(:,k),DOF_mat1(:,k),DOF_mat3(:,k)];
%     edge_nodes{k} = mesh_data.edges(dof_data.dofs_to_edges(DOF_mat{k})
end

temp1        = nonzeros(DOF_mat1); % Create temporary column vector
temp2        = nonzeros(DOF_mat2);
temp3        = nonzeros(DOF_mat3);

% Create vector containing nodes of each edge:
edge_nodes_1 = mesh_data.edges(dof_data.dofs_to_edges(temp1(1:2:end,1)),:); % Edges on contour
edge_nodes_2 = mesh_data.edges(dof_data.dofs_to_edges(temp2(1:2:end,1)),:); % First diagonal
edge_nodes_3 = mesh_data.edges(dof_data.dofs_to_edges(temp3(1:2:end,1)),:); % Second diagonal

clear temp1 temp2 temp3

edge_vecs_1  = mesh_data.node_coords(edge_nodes_1(:,1),:) - mesh_data.node_coords(edge_nodes_1(:,2),:); % Vector of edges
edge_vecs_2  = mesh_data.node_coords(edge_nodes_2(:,1),:) - mesh_data.node_coords(edge_nodes_2(:,2),:);
edge_vecs_3  = mesh_data.node_coords(edge_nodes_3(:,1),:) - mesh_data.node_coords(edge_nodes_3(:,2),:);
% PlotTriangleMeshRaw(mesh_data.node_coords,edge_nodes_3,1);

lim1 = length(edge_vecs_1);
lim2 = length(edge_vecs_2);
lim3 = length(edge_vecs_3);

% Determine which columns correspond, so that the correct angle can be
% determined.
% This is so that the MBF_mat values can be projected to the normal of the
% edge
[~,col1] = find(DOF_mat1);
[~,col2] = find(DOF_mat2);
lim1_for_2 = ismember(col1,col2);
lim1_for_2 = lim1_for_2(1:2:end);
[~,col2] = find(DOF_mat3);
lim1_for_3 = ismember(col1,col2);
lim1_for_3 = lim1_for_3(1:2:end);

theta_1      =  abs(90 - acosd(dot(edge_vecs_1(lim1_for_2,:),edge_vecs_2(1:lim2,:),2)./(vecnorm(edge_vecs_1(lim1_for_2,:),2,2).*vecnorm(edge_vecs_2(1:lim2,:),2,2))));
theta_2      =  abs(90 - acosd(dot(edge_vecs_1(lim1_for_3,:),edge_vecs_3(1:lim3,:),2)./(vecnorm(edge_vecs_1(lim1_for_3,:),2,2).*vecnorm(edge_vecs_3(1:lim3,:),2,2))));


B2_first_node = zeros(length(edge_nodes_2),1);
B3_second_node = zeros(length(edge_nodes_3),1);
if cyl_def.firstNode == "endCap" % The first endcap's nodes sit on the wrong side for some reason (edge_nodes defined like that)
     edge_nodes_2(1:numVertices,[1 2]) =   edge_nodes_2(1:numVertices,[2 1]);
end


for MBF_num = 1:3 % unity,sine,cosine
      maxVal = 0;
    if MBF_num > 1  && (oneEndcap || twoEndcaps)
        temp1 = MBF_mat(1:numVertices,MBF_num);  % Get max node of sine/cosine
%         if MBF_num == 2
             temp2 = MBF_mat(length(edge_nodes_2)-numVertices+1:length(edge_nodes_2),MBF_num); % Last endcap
%         else
%             temp2 = MBF_mat(length(edge_nodes_3)-numVertices+1:length(edge_nodes_3),MBF_num); % Last endcap
%         end
        [~, maxNode1] = max(temp1,[],1, 'linear');
        [maxVal, maxNode2] = max(temp2,[],1, 'linear');
        
        refVec2 = edge_vecs_1(maxNode1,:);
        refVec3 = edge_vecs_1(maxNode2,:);
        for k = 1:numVertices
            if twoEndcaps % The lack of (abs) here is due to the minus sign in the formula (cos(phi)ap - sin(phi)aphi)
                 theta_1(k) =  (90 - acosd(dot(refVec2,edge_vecs_2(k,:),2)./(vecnorm(refVec2,2,2).*vecnorm(edge_vecs_2(k,:),2,2))));
                 theta_2(end-numVertices+k) =(90 - acosd(dot(refVec3,edge_vecs_3(end-numVertices+k,:),2)./(vecnorm(refVec3,2,2).*vecnorm(edge_vecs_3(end-numVertices+k,:),2,2))));
            elseif cyl_def.firstNode == "endCap"
                 theta_1(k) =  (90 - acosd(dot(refVec2,edge_vecs_2(k,:),2)./(vecnorm(refVec2,2,2).*vecnorm(edge_vecs_2(k,:),2,2))));
            elseif cyl_def.lastNode == "endCap"
                 theta_2(end-numVertices+k) =(90 - acosd(dot(refVec3,edge_vecs_3(end-numVertices+k,:),2)./(vecnorm(refVec3,2,2).*vecnorm(edge_vecs_3(end-numVertices+k,:),2,2))));
            end
        end
    end

    if cyl_def.firstNode == "endCap" % This just makes the B2 matrix have the correct MBF value at the centre vertex
        B2_first_node = [repelem(maxVal,numVertices)';zeros(length(edge_nodes_2)-numVertices,1)];
%         B2_first_node = [MBF_mat(edge_nodes_2(1:numVertices,1),MBF_num);zeros(length(edge_nodes_2)-numVertices,1)];
    end
    if cyl_def.lastNode == "endCap"
        B3_second_node = [zeros(length(edge_nodes_3)-numVertices,1);repelem(maxVal,numVertices)'];
%          B3_second_node = [zeros(length(edge_nodes_3)-numVertices,1);MBF_mat(edge_nodes_3(end-numVertices+1:end,2),MBF_num)];
    end
    
    B1(1:2,:) = [MBF_mat(edge_nodes_1(:,1),MBF_num),MBF_mat(edge_nodes_1(:,2),MBF_num)]';
    B2(1:2,:) = [B2_first_node,MBF_mat(edge_nodes_2(:,2),MBF_num)]';
    if cyl_def.firstNode == "endCap"
        B2(2,1:numVertices) = B2(1,1:numVertices); % This just for endcap case
    end
    B3(1:2,:) = [MBF_mat(edge_nodes_3(:,1),MBF_num),B3_second_node]';
    if cyl_def.lastNode == "endCap"
        B3(1,end-numVertices+1:end) = B3(2,end-numVertices+1:end);
    end
    if MBF_num == 1
        % Don't want to include edges on the endcap for the unity case
        % (potloodpunt)
        if cyl_def.firstNode == "endCap"
            B2(1:2,1:numVertices) = [zeros(numVertices,1),zeros(numVertices,1)]';
        end
        if cyl_def.lastNode == "endCap"
            B3(1:2,end+1-numVertices:end) = [zeros(numVertices,1),zeros(numVertices,1)]';
        end
    end
    % Assume B1 (straight edges) are always aligned with axis (edge normal
    % is on axis)
    B2(:,:) = B2(:,:) .* [sind(theta_1)';sind(theta_1)']; % Get component of MBF_mat on edge normal
    B3(:,:) = B3(:,:) .* [sind(theta_2)';sind(theta_2)'];
    
% Solve 2x2 linear system
    for i = 1:lim1
        X1(:,i) = Rho(:,:,i)\B1(:,i);
    end
    for i = 1:lim2 % The diagonal sides have potentially one MBF less if there are endcaps/connections
        X2(:,i) = (Rho2(:,:,i))\B2(:,i);
    end
    for i = 1:lim3
        X3(:,i) = (Rho3(:,:,i))\B3(:,i);
    end
    
    col_iter = 1;
    node2 = 1;
    node3 = 1;
    xdom1 = 0;
    for MBF_node = 1:numNodes_new

        col_index = col_iter + (MBF_num-1);
        col_iter = col_iter + numMBF;

        % Skip the zero columns
        if DOF_mat1(1,MBF_node) ~= 0
            xdom_prev1 = xdom1 + 1;
            xdom1 = xdom1 + numVertices; % Next set of values
            U_Mat(DOF_mat1(1:2:end,MBF_node),col_index) = X1(1,xdom_prev1:xdom1); % RWG
            U_Mat(DOF_mat1(2:2:end,MBF_node),col_index) = X1(2,xdom_prev1:xdom1); % Linear
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