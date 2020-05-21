function [node_coords, elements] = ConnectToPlate(polygonPoints, node_coords_cyl, elements_cyl, numVertices, cyl_def)

iter = 0;
coord_iter = -2;

if cyl_def.firstNode == "conn"
    
   firstNodeCoords = node_coords_cyl(1:numVertices,:);
    
   num_poly = length(polygonPoints);
   num_cyl  = length(firstNodeCoords);
    
    
    if num_cyl > num_poly
        first_iter = num_poly;
        num_conn = num_cyl/num_poly; % Number of extra connection points
        firstCase = true;
    else
       first_iter = num_cyl; 
       num_conn = num_poly/num_cyl; % Number of extra connection points
       firstCase = false;
    end
    
    
    for i = 1:first_iter % Iterate smaller set
        for j = 0:num_conn-1 % Every point in smaller set has num_conn connections to larger set
            iter = iter + 1;
            coord_iter = coord_iter + 3;
            % Connect the points
            if firstCase
                connection_vecs(iter,:) = polygonPoints(i,:) - firstNodeCoords(i+j,:);
                % The problem is, we don't know how the points align, and
                % we can assume polygonPoints(i) is not 'under' or aligned
                % with firstNodeCoords(i)
                node_coords(coord_iter:coord_iter+2,:) = [polygonPoints(i,:); firstNodeCoords(i+j,:); firstNodeCoords(mod(i+j-1,numVertices)+1,:)];
%                 node_coords(coord_iter+3:coord_iter+5,:) = [polygonPoints(i,:); polygonPoints(i+j,:); firstNodeCoords(mod(i+j-1,numVertices)+1,:)];
                elements(iter, :) = [coord_iter, coord_iter+1, coord_iter+2, coord_iter+2];
            else
                connection_vecs(iter,:) = firstNodeCoords(i,:) - polygonPoints(i+j,:);
                node_coords(coord_iter:coord_iter+2,:) = [firstNodeCoords(i,:); polygonPoints(i+j,:); polygonPoints(mod(i+j-1,numVertices)+1,:)];
%                 node_coords(coord_iter+3:coord_iter+5,:) = [polygonPoints(i,:); polygonPoints(i+j,:); firstNodeCoords(mod(i+j-1,numVertices)+1,:)];
                elements(iter,:) = [coord_iter, coord_iter+1, coord_iter+2, coord_iter+2];
            end
        end
    end
    
end

if cyl_def.lastNode == "conn"
    
    lastNodeCoords = node_coords_cyl(end-numVertices+1:end,:);
    
    num_poly = length(polygonPoints);
    num_cyl  = length(lastNodeCoords);
    
    
    if num_cyl > num_poly
        first_iter = num_poly;
        num_conn = num_cyl/num_poly; % Number of extra connection points
        firstCase = true;
    else
       first_iter = num_cyl; 
       num_conn = num_poly/num_cyl; % Number of extra connection points
       firstCase = false;
    end
    
    
    for i = 1:first_iter % Iterate smaller set
        for j = 0:num_conn-1 % Every point in smaller set has num_conn connections to larger set
            % Connect the points
            if firstCase
                connection_vecs(i) = polygonPoints(i) - lastNodeCoords(i+j);
            else
                connection_vecs(i) = lastNodeCoords(i) - polygonPoints(i+j);
            end
        end
    end
    elements = connection_vecs;
    node_coords = connection_vecs;
    
    
    
    
end




end