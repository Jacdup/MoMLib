function [U_Mat, DOF_mat1, DOF_mat2, DOF_mat3] = MBF_Axial(mesh_data, dof_data, numVertices ,numMBF, numNodes, triangle_blah, cyl_def)

% The MBF is agnostic to signs on the basis functions, and only cares about
% if all the DOFS where the MBF is defined have the same orientation

% -------------------------------------------------------------------------
% Init
% -------------------------------------------------------------------------
oneEndcap = (cyl_def.firstNode == "endCap" && cyl_def.lastNode ~= "endCap") || (cyl_def.firstNode ~= "endCap" && cyl_def.lastNode == "endCap"); % Only one is an endcap
twoEndcaps = cyl_def.firstNode == "endCap" && cyl_def.lastNode == "endCap";
connection = cyl_def.firstNode == "conn" && cyl_def.lastNode == "conn";
    
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
    if cyl_def.coupling
        numNodes_new = numNodes + 4;% A total of 4 connections
    end
else
    numNodes_new = numNodes;
    endCapExclude = 0;
    connCapExclude = 0;
end

numDofs = size(dof_data.basis_supports,1);
numMBFNodes = numVertices*(numNodes_new) ;

U_Mat = zeros(numDofs, numNodes_new*numMBF);
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
%     if cyl_def.coupling
%         contour_nodes = (1:(numNodes+4)*numVertices)';
%     end
    contour_nodes = [contour_nodes; cyl_def.last_element_val +  cyl_def.plate_polygon_nodes;  cyl_def.last_element_val + cyl_def.plate_polygon_nodes_end];
elseif cyl_def.firstNode == "conn"
    contour_nodes(end-numVertices+1:end) =   ((numNodes+1)*numVertices) + 1 + cyl_def.plate_polygon_nodes ; % Since the polygon nodes are not sequential anymore
elseif cyl_def.lastNode == "conn"
    contour_nodes = [contour_nodes; cyl_def.last_element_val + cyl_def.plate_polygon_nodes_end];
end

% Contour_nodes can now be used to set the values through linear indexing:
MBF_mat = zeros(numMBFNodes*2, 3); % Don't really care how big this matrix is, as long as the rows correspond to the nodes.
MBF_mat(contour_nodes, :) = [ones_mat,sin_mat',cos_mat'];


% MBF_mat = MBF_mat(1:end-36,:); % Temporary for coupling

row = -1;
iter = 1;
col = 1;
linear_row = 0;
rho_ind = 0;



for numCyl = 1:cyl_def.coupling+1 % Iterates twice if there are two cylinders
    
    DOF_mat1 = zeros(numVertices*2,numNodes_new);
    DOF_mat2 = zeros(numVertices*2,numNodes_new);
    DOF_mat3 = zeros(numVertices*2,numNodes_new);
    X1       = zeros(2,numMBFNodes);
    X2       = zeros(2,numMBFNodes);
    X3       = zeros(2,numMBFNodes);
    
    
%     len_tri_mat = 1:2:(length(triangle_blah)-vert_num-1-endCapExclude-cyl_def.num_plate_nodes - connCapExclude); % The length of the matrix with just the cylinder(s)
    startNode = 1;
    endNode = (length(triangle_blah)-vert_num-1-endCapExclude-cyl_def.num_plate_nodes - connCapExclude);
     
%     if cyl_def.coupling
%         startNode = 1;
%         endNode = 1399;
%         if numCyl == 2
% %             startNode = endNode + 1;
%             startNode = 1401;
%             endNode = 2839;
%         end
% %         endNode = round((length(triangle_blah)-vert_num-1-endCapExclude-cyl_def.num_plate_nodes - connCapExclude)/(mod(numCyl,2)+1));
% %         endNode = (length(triangle_blah) -1- cyl_def.num_plate_nodes)/(mod(numCyl,2)+1);
%     end
    len_tri_mat = startNode:2:endNode;
    
    Rho = [1,1;1,-1]; % This is the matrix of RWG and linear components at the two edge nodes
    Rho = repmat(Rho,1,1,numMBFNodes);
    Rho2 = Rho;
    Rho3 = Rho;

    
    % Select DOFs associated with the MBF
    for i = len_tri_mat
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

%         if i <= numVertices*2 && cyl_def.firstNode == "conn" % Assign first MBF with new DOFs from generated connection
        if (i <= startNode + vert_num) && cyl_def.firstNode == "conn" % Assign first MBF with new DOFs from generated connection
            % This is NOT the transition from cyl -> plate
            DOF_mat1(row:row+1,col) = [triangle_blah(i+1,9);triangle_blah(i+1,15)];
            DOF_mat2(row:row+1,col) = [triangle_blah(i+1-last,8);triangle_blah(i+1-last,14)];
        end
        if (i >= endNode - vert_num +1) && cyl_def.lastNode == "conn"  % Assign last MBF with new DOFs from generated connection
            % This is NOT the transition from cyl -> plate
            DOF_mat3(row:row+1,col) = [triangle_blah(i+vert_num+2-last,8);triangle_blah(i+vert_num+2-last,14)];
        end
        
        if last == 1
            row = -1;
            linear_row = 0;
            col = col + 1;
            if col == 36
                test = 1;
            end
        end
    end
    for i = numVertices:numVertices:length(Rho) % DOFs at the end of each node have switched signs
        Rho2(:,:,i) = Rho2(:,:,i) .* [1,-1;1,-1];
        Rho3(:,:,i) = Rho3(:,:,i) .* [1,-1;1,-1];
        %     Rho(:,:,i) = Rho(:,:,i) .* [1,-1;1,-1];
    end

    % -------------------------------------------------------------------------
    % -----------------------------Endcap stuff--------------------------------
    % -------------------------------------------------------------------------
    if oneEndcap || twoEndcaps || connection
        row1 = -1; % Fill up a new column, from row 1
        row2 = -1;
        linear_row = 0;
        if cyl_def.firstNode == "endCap" || cyl_def.firstNode == "conn"
%             for i = 1:2:(vert_num) % First endcap, or connection triangles
            for i = startNode:2:startNode + vert_num
                linear_row = linear_row + 1;
                row1 = row1+ 2;
                last = 0;
                if (i == vert_num)
                    last = 1;
                end
                
                if cyl_def.firstNode == "conn"
                    DOF_mat1(row1:row1+1, col) = [triangle_blah(i,7); triangle_blah(i,13)];
                    DOF_mat3(row1:row1+1, col) = [triangle_blah(i+1,8-last); triangle_blah(i+1,14-last)];
                    
                    %                  Rho2(:,:, Rho_index) = Rho2(:,:,Rho_index) .* [-1,-1;-1,-1];
                    if connection
                        Rho_index = sub2ind(size(DOF_mat1(1:2:end,:)), linear_row, col);
                        Rho3(:,:, Rho_index) = Rho3(:,:,Rho_index) .* [-1,-1;-1,-1]; % This fixes the 'connection' case (2 connections)
                    end
                    %
                else
                    DOF_mat1(row1:row1+1, col) = [triangle_blah(i,9); triangle_blah(i,15)];
                    DOF_mat3(row1:row1+1, col) = [triangle_blah(i,7+last); triangle_blah(i,13+last)];
                end
                
            end
            %         DOF_mat1(:,col) = sort(DOF_mat1(:,col));
        end
        
        linear_row = 0;
        if cyl_def.lastNode == "endCap" || cyl_def.lastNode == "conn"
%             for i = (length(triangle_blah)-vert_num-endCapExclude- connCapExclude - cyl_def.num_plate_nodes+1):2:length(triangle_blah)-endCapExclude-connCapExclude-cyl_def.num_plate_nodes % Second endcap, or connection triangles
            for i = endNode+2:2:endNode+vert_num+1
                row2 = row2 + 2;
                linear_row = linear_row + 1;
                last = 0;
                
%                 if i ==  length(triangle_blah)-endCapExclude-connCapExclude-cyl_def.num_plate_nodes
                if i == endNode+vert_num+1
                    last = 1;
                end
                Rho_index = sub2ind(size(DOF_mat1(1:2:end,:)), linear_row, col+1+extra_dof_col);
                %             DOF_mat1(row2:row2+1, col+1+extra_dof_col) = [triangle_blah(i,7); triangle_blah(i,13)];
                %             if cyl_def.lastNode == "conn"
                if cyl_def.firstNode == "conn" && cyl_def.lastNode ~= "conn" % Dofs are assigned differently when first node is a connection
                    DOF_mat2(row2:row2+1, col+1+extra_dof_col) = [triangle_blah(i+last,7+last); triangle_blah(i+last,13+last)];
                    DOF_mat1(row2:row2+1, col+1+extra_dof_col) = [triangle_blah(i+1,7); triangle_blah(i+1,13)];
                elseif cyl_def.lastNode == "conn"
                    DOF_mat1(row2:row2+1, col+1+extra_dof_col) = [triangle_blah(i,7); triangle_blah(i,13)];
                    DOF_mat2(row2:row2+1, col+1+extra_dof_col) = [triangle_blah(i+1-last,8); triangle_blah(i+1-last,14)];
                    if connection
                        Rho_index = sub2ind(size(DOF_mat1(1:2:end,:)), linear_row, col+1+extra_dof_col);
                        Rho(:,:, Rho_index) = Rho(:,:,Rho_index) .* [-1,-1;-1,-1]; % This fixes the 'connection' case (2 connections)
                        %                     Rho2(:,:, Rho_index) = Rho2(:,:,Rho_index) .* [-1,-1;-1,-1];
                        %                     Rho3(:,:, Rho_index) = Rho3(:,:,Rho_index) .* [-1,-1;-1,-1];% uncomment this for pretty current
                    end
                    %                 Rho3(:,:, Rho_index) = Rho3(:,:,Rho_index) .* [-1,-1;-1,-1];
                    %                 Rho2(:,:, Rho_index) = Rho2(:,:,Rho_index) .* [1,-1;1,-1];
                else
                    DOF_mat1(row2:row2+1, col+1+extra_dof_col) = [triangle_blah(i,7); triangle_blah(i,13)];
                    DOF_mat2(row2:row2+1, col+1+extra_dof_col) = [triangle_blah(i,9-last); triangle_blah(i,15-last)];
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
    
    if (cyl_def.firstNode == "endCap" && cyl_def.lastNode == "endCap") || (cyl_def.firstNode == "conn" && cyl_def.lastNode ~= "conn")
        DOF_mat1 = [circshift(DOF_mat1(:,1:end-1), [0 1]), DOF_mat1(:,end)]; % Swap columns, so that DOFs are ascending from column 1
        DOF_mat2 = [circshift(DOF_mat2(:,1:end-1), [0 1]), DOF_mat2(:,end)];
        DOF_mat3 = [circshift(DOF_mat3(:,1:end-1), [0 1]), DOF_mat3(:,end)];
        %     DOF_mat1(:,end) = sort(DOF_mat1(:,end));
        %     DOF_mat2(:,end) = sort(DOF_mat2(:,end));
        %     DOF_mat3(:,end) = sort(DOF_mat3(:,end));
        Rho(:,:,1:numVertices) = Rho(:,:,1:numVertices) .* [-1,-1;-1,-1]; % I REALLY don't know why this is suddenly necessary (since 25/05/2020)
        %      Rho2(:,:,1:numVertices) = Rho2(:,:,1:numVertices) .* [-1,-1;-1,-1]; % I REALLY don't know why this is suddenly necessary (since 25/05/2020)
    elseif cyl_def.firstNode == "endCap"
        Rho(:,:,end-numVertices+1:end) = Rho(:,:,end-numVertices+1:end) .* [-1,-1;-1,-1]; % I REALLY don't know why this is suddenly necessary (since 25/05/2020)
    end
    
    
    temp1        = nonzeros(DOF_mat1); % Create temporary column vector
    temp2        = nonzeros(DOF_mat2);
    temp3        = nonzeros(DOF_mat3);
    % Create vector containing nodes of each edge:
    edge_nodes_1 = mesh_data.edges(dof_data.dofs_to_edges(temp1(1:2:end,1)),:); % Edges on contour
    edge_nodes_2 = mesh_data.edges(dof_data.dofs_to_edges(temp2(1:2:end,1)),:); % First diagonal
    edge_nodes_3 = mesh_data.edges(dof_data.dofs_to_edges(temp3(1:2:end,1)),:); % Second diagonal
    B1       = zeros(2,length(edge_nodes_1(:,1)));
    B2       = zeros(2,length(edge_nodes_2(:,2)));
    B3       = zeros(2,length(edge_nodes_3(:,1)));
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
    [~,col1] = find(DOF_mat1);
    
    [~,col2] = find(DOF_mat2);
    lim1_for_2 = ismember(col1,col2);
    lim1_for_2 = lim1_for_2(1:2:end);
    
    [~,col2] = find(DOF_mat3);
    lim1_for_3 = ismember(col1,col2);
    lim1_for_3 = lim1_for_3(1:2:end);
    %
    % lim1_for_3 = 1:lim3;
    % lim1_for_2 = 1:lim2;
    
    
%     theta_1      =  abs(90 - acosd(dot(edge_vecs_1(lim1_for_2,:),edge_vecs_2(1:lim2,:),2)./(vecnorm(edge_vecs_1(lim1_for_2,:),2,2).*vecnorm(edge_vecs_2(1:lim2,:),2,2))));
%     theta_2      =  abs(90 - acosd(dot(edge_vecs_1(lim1_for_3,:),edge_vecs_3(1:lim3,:),2)./(vecnorm(edge_vecs_1(lim1_for_3,:),2,2).*vecnorm(edge_vecs_3(1:lim3,:),2,2))));
    
    theta_1      =  abs(acosd(dot(edge_vecs_1(lim1_for_2,:),edge_vecs_2(1:lim2,:),2)./(vecnorm(edge_vecs_1(lim1_for_2,:),2,2).*vecnorm(edge_vecs_2(1:lim2,:),2,2))));
    theta_2      =  abs(acosd(dot(edge_vecs_1(lim1_for_3,:),edge_vecs_3(1:lim3,:),2)./(vecnorm(edge_vecs_1(lim1_for_3,:),2,2).*vecnorm(edge_vecs_3(1:lim3,:),2,2))));
    for MBF_num = 1:3
        
        B1(1:2,:) = [MBF_mat(edge_nodes_1(:,1),MBF_num),MBF_mat(edge_nodes_1(:,2),MBF_num)]';
        B2(1:2,:) = [zeros(lim2,1),MBF_mat(edge_nodes_2(:,2),MBF_num)]';
        B3(1:2,:) = [MBF_mat(edge_nodes_3(:,1),MBF_num),zeros(lim3,1)]';
%         B2 = B2(:,:) .* [sind(theta_1)';sind(theta_1)'];
%         B3 = B3(:,:) .* [sind(theta_2)';sind(theta_2)'];
        B2 = B2(:,:) .* [cosd(2.*theta_1)';cosd(2.*theta_1)'];
        B3 = B3(:,:) .* [cosd(2.*theta_2)';cosd(2.*theta_2)'];
        
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
end