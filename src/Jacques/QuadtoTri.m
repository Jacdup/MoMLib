function [new_triangles,triangles] = QuadtoTri(element)

num_nodes = length(element);

triangle1 = zeros(num_nodes,3);
triangle2 = zeros(num_nodes,3);

for node = 1:num_nodes
    % Define triangle corners
    c1 = element(node,1);
    c2 = element(node,2);
    c3 = element(node,4);
    c4 = element(node,3);   
    triangle1(node,1:3) = [c1 c2 c3];
    triangle2(node,1:3) = [c2 c3 c4];
end

triangles = zeros(length(triangle1)*2,3);
triangles(1:2:end-1,1:3) = triangle1;
triangles(2:2:end,1:3)   = triangle2;
new_triangles = triangles;
% We have to convert to format used by Robey:
% new_triangles = {};
% for node = 1:num_nodes*2
%     new_triangles = [new_triangles;triangles(node,1:3),[0 0 0 0 0 0]]; 
% end
