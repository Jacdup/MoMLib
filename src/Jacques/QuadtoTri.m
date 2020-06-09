function [new_triangles,triangles] = QuadtoTri(element, numVertices, cyl_def)

num_nodes = length(element);
num_diff = 0;
connection_flag = 0;
both = (cyl_def.firstNode == "endCap" && cyl_def.lastNode == "endCap") || (cyl_def.firstNode == "conn" && cyl_def.lastNode == "conn");
oneEach = (cyl_def.firstNode == "endCap" && cyl_def.lastNode == "conn") || (cyl_def.firstNode == "conn" && cyl_def.lastNode == "endCap");
one = cyl_def.firstNode == "endCap" || cyl_def.lastNode == "endCap";
% if connection_flag
%    num_diff = (6*numVertices); 
% end
% if connection_flag == 1 && endcap == 0
%    num_diff = 2 * numVertices; 
% end
if both
    num_diff = 2*numVertices;
%     num_diff = 4*numVertices; % temporary
elseif oneEach
%     num_diff = 3*numVertices; % Previous algo 
   num_diff = numVertices; 
elseif one
    num_diff = numVertices;
else
   num_diff = 0; 
end
num_diff = 0;

% if cyl_def.firstNode == "conn" || cyl_def.lastNode == "conn"
%         num_nodes = num_nodes - num_diff; % Since the last elements are already triangles
%     triangle3 = zeros(num_diff,3);
%     
%     for node = num_nodes+1:length(element)
%         c1 = element(node,1);
%         c2 = element(node,2);
%         c3 = element(node,4);
%         %         c4 = element(node,3);
%         triangle3(node-num_nodes,1:3) = [c1 c2 c3];
%         
%     end
% end

% if connection_flag
%     num_diff = (2*numVertices);
%     num_nodes = num_nodes - num_diff; % Since the last elements are already triangles
%     
%     for node = num_nodes+1:length(element)
%         c1 = element(node,1);
%         c2 = element(node,2);
%         c3 = element(node,4);
%         %         c4 = element(node,3);
%         triangle3(node-num_nodes,1:3) = [c1 c2 c3];
%     end
% end

% if oneEach || both || one
% %     num_diff = (2*numVertices);
%     num_nodes = num_nodes - num_diff; % Since the last elements are already triangles
%     triangle3 = zeros(num_diff,3);
%     
%     for node = num_nodes+1:length(element)
%         c1 = element(node,1);
%         c2 = element(node,2);
%         c3 = element(node,4);
%         %         c4 = element(node,3);
%         triangle3(node-num_nodes,1:3) = [c1 c2 c3];
%         
%     end
% end
triangle1 = zeros(num_nodes,3);
triangle2 = zeros(num_nodes,3);

for node = 1:num_nodes
    % Define triangle corners
    c1 = element(node,1);
    c2 = element(node,2);
    c3 = element(node,4);
    c4 = element(node,3);
    triangle1(node,1:3) = [c1 c2 c3];
    if c3 ~= c4 % Check for degeneracy
        triangle2(node,1:3) = [c2 c3 c4];
    end
%     if mod(node,numVertices-1) == 0
%         triangle1(node,1:3) = [c2 c1 c3];
%         triangle2(node,1:3) = [c3 c2 c4];
%     end
end


triangles = zeros(length(triangle1)*2 + num_diff,3);
triangles(1:2:end-1-num_diff,1:3) = triangle1;
triangles(2:2:end-num_diff,1:3)   = triangle2;
% if cyl_def.lastNode == "endCap" || cyl_def.firstNode == "endCap"
%     triangles(end-num_diff+1:end, 1:3) = triangle3;
% end
new_triangles = triangles;

new_triangles( ~any(new_triangles,2), : ) = [];  %rows
% We have to convert to format used by Robey:
% new_triangles = {};
% for node = 1:num_nodes*2
%     new_triangles = [new_triangles;triangles(node,1:3),[0 0 0 0 0 0]];
% end
