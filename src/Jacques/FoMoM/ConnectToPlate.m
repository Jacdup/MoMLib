function [node_coords, elements] = ConnectToPlate(polygonPoints,plate_polygon_nodes, node_coords_cyl, elements_cyl,numVertices, cyl_def)

iter = 0;
elem_iter = -1;
coord_iter = 0;
last_element_val = max(max(elements_cyl));
node_coords = polygonPoints;

if cyl_def.firstNode == "conn"
    
   firstNodeCoords = node_coords_cyl(1:numVertices,:);
    
   num_poly = length(polygonPoints); % Number of polygon vertices
   num_cyl  = length(firstNodeCoords); % Number cylinder vertices
    
    
    if num_cyl > num_poly
        first_iter = num_poly;
        num_conn = num_cyl/num_poly; % Number of extra connection points
        firstCase = true;
    else
       first_iter = num_cyl; 
       num_conn = num_poly/num_cyl; % Number of extra connection points
       firstCase = false;
    end
    max_coord_iter = num_conn * first_iter;
    % The problem is, we don't know how the points align, and
    % we can assume polygonPoints(i) is not 'under' or aligned
    % with firstNodeCoords(i)
    % SO, first find points that align
    min = 1000;
    align_offset = 1;
    for i = 1:first_iter % Cycle through points and find minimum distance. This will be the index that aligns the cylinder with the polygon.
        if firstCase
             vec_len = firstNodeCoords(1,:) - polygonPoints(i,:);
        else
            vec_len =  polygonPoints(1,:) - firstNodeCoords(i,:) ;
        end
        if norm(vec_len) < min
            min = norm(vec_len);
            align_offset = 0;
        end
    end
    
    for i = 1:first_iter % Iterate smaller set
        for j = 0:num_conn-1 % Every point in smaller set has num_conn connections to larger set
            iter = iter + 1;
            elem_iter = elem_iter + 2;
            coord_iter = coord_iter + 1;
            % Connect the points
            if firstCase % Cylinder has more points
                
                node_coords(coord_iter:coord_iter+2,:) = [polygonPoints(mod(num_poly - i,num_poly)+1,:); firstNodeCoords(mod(i-j+align_offset,num_cyl)+1,:); firstNodeCoords(mod(i-j-1+align_offset,num_cyl)+1,:)];
                %                 node_coords(coord_iter+3:coord_iter+5,:) = [polygonPoints(i,:); polygonPoints(i+j,:); firstNodeCoords(mod(i+j-1,numVertices)+1,:)];
                elements(elem_iter, :) = [coord_iter, coord_iter+1, coord_iter+2, coord_iter+2];
                elements(elem_iter + 1, :) = [coord_iter+1, coord_iter+2,  mod(coord_iter+3, max_coord_iter)+1,  mod(coord_iter+3, max_coord_iter)+1];
            else % PolygonPoints has more points (or they are equal)
                %                 elements(elem_iter,:) = [elem_iter +last_element_val, elem_iter +1 + last_element_val,
                %                 node_coords(coord_iter:coord_iter+1,:) = [polygonPoints(mod(i-j+align_offset,num_poly)+1,:); polygonPoints(mod(i-j-1+align_offset,num_poly)+1,:)];
                %                 node_coords(coord_iter+3:coord_iter+5,:) = [polygonPoints(i,:); polygonPoints(i+j,:); firstNodeCoords(mod(i+j-1,numVertices)+1,:)];
                %                 elements(elem_iter,:) = [last_element_val + plate_polygon_nodes(coord_iter) ,last_element_val + plate_polygon_nodes(mod(coord_iter,length(plate_polygon_nodes))+1), mod(num_cyl - coord_iter + align_offset,num_cyl)+1, mod(num_cyl - coord_iter+align_offset,num_cyl)+1];
                %                 elements(elem_iter+1,:) = [mod(num_cyl - coord_iter +align_offset,num_cyl)+1, mod(num_cyl - coord_iter + 1+ align_offset, num_cyl)+1, plate_polygon_nodes(coord_iter)+last_element_val, plate_polygon_nodes(coord_iter)+last_element_val];
                T1_cnr1 = last_element_val + plate_polygon_nodes(coord_iter);
                T1_cnr2 = last_element_val + plate_polygon_nodes(mod(coord_iter,length(plate_polygon_nodes))+1);
%                 T1_cnr3 = mod(coord_iter + align_offset,num_cyl)+1;
                
                T2_cnr1 = mod(coord_iter +align_offset,num_cyl)+1;
                T2_cnr2 = mod(coord_iter - 1+ align_offset, num_cyl)+1;
%                 T2_cnr3 = plate_polygon_nodes(coord_iter)+last_element_val;
                
                elements(elem_iter  ,:) = [T1_cnr1, T1_cnr2, T2_cnr2, T2_cnr2]; % Triangle 1 (Meshing degenerate quads actually)
                elements(elem_iter+1,:) = [T2_cnr1, T2_cnr2, T1_cnr2, T1_cnr2]; % Triangle 2
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