function [DOF_mat1,DOF_mat2,DOF_mat3] = getAxialDOFs(tri_mat,numVertices,num_interior_nodes, len_tri_mat)
DOF_mat1 = zeros(numVertices*2,num_interior_nodes+2);
DOF_mat2 = zeros(numVertices*2,num_interior_nodes+2);
DOF_mat3 = zeros(numVertices*2,num_interior_nodes+2);
vert_num = (2*numVertices)-1;

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
        DOF_mat1(row:row+1,col) = [tri_mat(i+1,7);tri_mat(i+1,13)];
        DOF_mat2(row:row+1,col) = [tri_mat(i+last,7+last);tri_mat(i+last,13+last)];
        DOF_mat3(row:row+1,col) = [tri_mat(i+vert_num+1,7+last);tri_mat(i+vert_num+1,13+last)];

        if i <= numVertices*2 && cyl_def.firstNode == "conn" % Assign first MBF with new DOFs from generated connection
            % This is NOT the transition from cyl -> plate
            DOF_mat1(row:row+1,col) = [tri_mat(i+1,9);tri_mat(i+1,15)];
            DOF_mat2(row:row+1,col) = [tri_mat(i+1-last,8);tri_mat(i+1-last,14)];
        end
        if cyl_def.lastNode == "conn" && i >= len_tri_mat - vert_num +1 % Assign last MBF with new DOFs from generated connection
            % This is NOT the transition from cyl -> plate
            DOF_mat3(row:row+1,col) = [tri_mat(i+vert_num+2-last,8);tri_mat(i+vert_num+2-last,14)];
        end

        if last == 1
            row = -1;
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
                    DOF_mat1(row1:row1+1, col) = [tri_mat(i,7); tri_mat(i,13)];
                    DOF_mat3(row1:row1+1, col) = [tri_mat(i+1,8-last); tri_mat(i+1,14-last)];

    %                  Rho2(:,:, Rho_index) = Rho2(:,:,Rho_index) .* [-1,-1;-1,-1];
                if connection
                    Rho_index = sub2ind(size(DOF_mat1(1:2:end,:)), linear_row, col);
                    Rho3(:,:, Rho_index) = Rho3(:,:,Rho_index) .* [-1,-1;-1,-1]; % This fixes the 'connection' case (2 connections)
                end
    %                  
                else
    %                  Rho_index = sub2ind(size(DOF_mat1(1:2:end,:)), linear_row, col);
    %                  Rho(:,:, Rho_index) = Rho(:,:,Rho_index) .* [-1,-1;-1,-1];
                    DOF_mat1(row1:row1+1, col) = [tri_mat(i,9); tri_mat(i,15)];
                    DOF_mat3(row1:row1+1, col) = [tri_mat(i,7+last); tri_mat(i,13+last)];
                end

            end
    %         DOF_mat1(:,col) = sort(DOF_mat1(:,col));
        end

        linear_row = 0;
        if cyl_def.lastNode == "endCap" || cyl_def.lastNode == "conn" 
            for i = (length(tri_mat)-vert_num-endCapExclude- connCapExclude - cyl_def.num_plate_nodes+1):2:length(tri_mat)-endCapExclude-connCapExclude-cyl_def.num_plate_nodes % Second endcap, or connection triangles
                row2 = row2 + 2;
                linear_row = linear_row + 1;
                last = 0;

                if i ==  length(tri_mat)-endCapExclude-connCapExclude-cyl_def.num_plate_nodes
                    last = 1;
                end
    %             Rho_index = sub2ind(size(DOF_mat1(1:2:end,:)), linear_row, col+1);
    %             DOF_mat1(row2:row2+1, col+1+extra_dof_col) = [tri_mat(i,7); tri_mat(i,13)];
    %             if cyl_def.lastNode == "conn" 
                if cyl_def.firstNode == "conn" && cyl_def.lastNode ~= "conn" % Dofs are assigned differently when first node is a connection
                   DOF_mat2(row2:row2+1, col+1) = [tri_mat(i+last,7+last); tri_mat(i+last,13+last)];
                    DOF_mat1(row2:row2+1, col+1) = [tri_mat(i+1,7); tri_mat(i+1,13)];
                elseif cyl_def.lastNode == "conn"
                    DOF_mat1(row2:row2+1, col+1) = [tri_mat(i,7); tri_mat(i,13)];
                    DOF_mat2(row2:row2+1, col+1) = [tri_mat(i+1-last,8); tri_mat(i+1-last,14)];
                    if connection
                        Rho_index = sub2ind(size(DOF_mat1(1:2:end,:)), linear_row, col+1);
                        Rho(:,:, Rho_index) = Rho(:,:,Rho_index) .* [-1,-1;-1,-1]; % This fixes the 'connection' case (2 connections)
    %                     Rho2(:,:, Rho_index) = Rho2(:,:,Rho_index) .* [-1,-1;-1,-1];
    %                     Rho3(:,:, Rho_index) = Rho3(:,:,Rho_index) .* [-1,-1;-1,-1];% uncomment this for pretty current
                    end
    %                 Rho3(:,:, Rho_index) = Rho3(:,:,Rho_index) .* [-1,-1;-1,-1];
    %                 Rho2(:,:, Rho_index) = Rho2(:,:,Rho_index) .* [1,-1;1,-1];
                else
                    DOF_mat1(row2:row2+1, col+1) = [tri_mat(i,7); tri_mat(i,13)];
                    DOF_mat2(row2:row2+1, col+1) = [tri_mat(i,9-last); tri_mat(i,15-last)];
                end



    %              Rho3(:,:, Rho_index) = Rho3(:,:,Rho_index) .* [-1,-1;-1,-1];
    %             [Rho(:,:,Rho_index), Rho2(:,:,Rho_index)] = getSigns_new(tri_mat,7,9-last,i,i);

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
         DOF_mat2 = [circshift(DOF_mat2(:,1:end-1), [0 1]), DOF_mat2(:,end)];
        row = -1;
        if cyl_def.firstNode == "endCap"
            for i = length(tri_mat)-endCapExclude+1-cyl_def.num_plate_nodes:1:length(tri_mat)-(numVertices*(cyl_def.lastNode == "endCap")-cyl_def.num_plate_nodes)% Every odd row, only first endcap
                row = row + 2;
                ind = ind + 1;
    %             Rho2(:,:,ind) = [-1,1;-1,-1];
    %            Add the edges on the end cap
                DOF_mat2(row:row+1,1) = [tri_mat(i,7);tri_mat(i,13)];
                if i == length(tri_mat)-(numVertices*(cyl_def.lastNode == "endCap")-cyl_def.num_plate_nodes)
                    Rho2(:,:,ind) = [-1,1;-1,-1];
                     DOF_mat2(row:row+1,1) = [tri_mat(i,8);tri_mat(i,14)];
                     col = col+1;
                end
            end
        end
    %      DOF_mat2(:,~any(DOF_mat,1)) = []; % Remove zero columns
    %     Rho(:,:,i+2) = [1,-1;1,1];
        row = -1;
        if cyl_def.lastNode == "endCap" && cyl_def.firstNode ~= "conn"
            for i = length(tri_mat)-cyl_def.num_plate_nodes:-1:length(tri_mat)-numVertices+1-cyl_def.num_plate_nodes% Second endcap
                row = row + 2; 
                if i == length(tri_mat)-numVertices+1-cyl_def.num_plate_nodes
                    Rho3(:,:,i) = [-1,-1;-1,1];
                    DOF_mat3(row:row+1,col) = [tri_mat(i,7);tri_mat(i,13)];%;tri_mat(i+1,8);tri_mat(i+1,14)]; 
                else
                    DOF_mat3(row:row+1,col) = [tri_mat(i,8);tri_mat(i,14)];  
                    if i == length(tri_mat)-cyl_def.num_plate_nodes
                         Rho3(:,:,i) = [1,-1;1,1];
                    end
                end

            end
        end
    end


    % -------------------------------------------------------------------------

    if (cyl_def.firstNode == "endCap" && cyl_def.lastNode == "endCap") || (cyl_def.firstNode == "conn" && cyl_def.lastNode ~= "conn")
        DOF_mat1 = [circshift(DOF_mat1(:,1:end-1), [0 1]), DOF_mat1(:,end)]; % Swap columns, so that DOFs are ascending from column 1
        DOF_mat3 = [circshift(DOF_mat3(:,1:end-1), [0 1]), DOF_mat3(:,end)];
        Rho(:,:,1:numVertices) = Rho(:,:,1:numVertices) .* [-1,-1;-1,-1]; % I REALLY don't know why this is suddenly necessary (since 25/05/2020)
    elseif cyl_def.firstNode == "endCap"
        Rho(:,:,end-numVertices+1:end) = Rho(:,:,end-numVertices+1:end) .* [-1,-1;-1,-1]; % I REALLY don't know why this is suddenly necessary (since 25/05/2020)
    end
    % Create cell array that has all the DOFs associated with all the MBFs 
    % (One MBF per cell)
    DOF_mat2( :, ~any(DOF_mat2,1) ) = [];  %remove zero columns
    DOF_mat3( :, ~any(DOF_mat3,1) ) = [];  %remove zero column
    for k = 1:size(DOF_mat1,2)
        DOF_mat{k} = [DOF_mat2(:,k),DOF_mat1(:,k),DOF_mat3(:,k)];
    %     edge_nodes{k} = mesh_data.edges(dof_data.dofs_to_edges(DOF_mat{k})
    end

end