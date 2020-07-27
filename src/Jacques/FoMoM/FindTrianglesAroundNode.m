function [tri_indices] = FindTrianglesAroundNode(contour_point, tri_nodes,tri_coords, num_vertices, coord, test_nodes)

% This function simply cycles through all the triangles and selects those
% who have two nodes on the MBF node (along the contour of a cylinder) to
% be excited.

% 2020-06-12. JT du Plessis

% EPS = 1e-6;
EPS = 0.005;
% tri_indices = zeros(num_vertices,2);

if ~isempty(test_nodes)
    iter = 0;
    % Get indices of tri_coords that are members of test_coords
    [~,LocB] = ismember(tri_coords,tri_coords(test_nodes,:), 'rows');
    
    % The triangle indices are the triangles that have two consecutive
    % test_coords
    for i = 1:2:length(test_nodes)
        
        edge = [test_nodes(i),test_nodes(i+1)];
        
        for j = 1:3
            % Find triangles that have this edge:
            [~,edge_tri_indices(:,j)] = ismember(tri_nodes(:,j),edge);
            [tri_indices(:,j),col(:,j)] = find(edge_tri_indices(:,j));
        end

    end
    
    
else % No connection, Simply check one coordinate value for finding indices around a node
    
    iter = 1;
    for i = 1:length(tri_nodes)
            
        LocB = abs(contour_point - tri_coords(tri_nodes(i,:),coord))<EPS; % Only look at 'coord'-coordinate here
        
        if length(nonzeros(LocB)) > 1 % If there are more than two points on the contour, then we know the triangle is part of the excitation/load
            tri_indices(iter,1) = i;
            
            if i == 1
                test = 1;
            end
            
            % Find corresponding (negative) tri that shares the edge:
            t = ismember(tri_nodes(:,:),tri_nodes(i,LocB));
            temp = find(sum(t~=0,2) > 1);
            tri_indices(iter,2) = temp(temp ~= i);
            
            
            iter = iter + 1;
        end
        
    end
    
    tri_indices = unique(sort(tri_indices,2),'rows'); % Keep only one copy of each triangle pair
    
end

end