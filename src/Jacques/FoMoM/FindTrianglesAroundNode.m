function [tri_indices] = FindTrianglesAroundNode(contour_point, tri_nodes,tri_coords)

% This function simply cycles through all the triangles and selects those
% who have two nodes on the MBF node (along the contour of a cylinder) to
% be excited.

% 2020-06-12. JT du Plessis

EPS = 1e-2;
EPS = 0.005;

iter = 1;
for i = 1:length(tri_nodes)
    
    LocB = abs(contour_point- tri_coords(tri_nodes(i,:),1))<EPS; % Only look at z-coordinate here
   
    if length(nonzeros(LocB)) > 1 % If there are more than two points on the contour, then we know the triangle is part of the excitation/load
       tri_indices(iter,1) = i; 
       
       % Find corresponding (negative) tri that shares the edge:
       t = ismember(tri_nodes(:,:),tri_nodes(i,LocB));
       temp = find(sum(t~=0,2) > 1);
       tri_indices(iter,2) = temp(temp ~= i);


       iter = iter + 1;
    end
    
end

tri_indices = unique(sort(tri_indices,2),'rows'); % Keep only one copy of each triangle pair


end