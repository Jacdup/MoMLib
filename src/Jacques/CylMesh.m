 function [points,quadElements, element] = CylMesh(Contour,num_vertices,rho, cyl_def)
% -------------------------------------------------------------------------
%     Create tube with planar quadrilateral meshing
%     V5.0 This function meshes a straight or curved section of wire
%     (defined with a 3D contour) with quadrilaterals. The mechanism for
%     rotating a point p about a unit vector r is quaternion multiplication
%     This function supports a maximum curve/angle of 90 degrees in the
%     cable
%
%       INPUT:
%       rho     ==>  Radius of shield
%       num_vertices  ==>  Number of vertices
%       Contour ==>  3D contour info (format: columns = x y z, rows =
%       nodes)
%
%       OUTPUT:
%       points  ==>  Point coordinates around nodes (format: columns = x y z, rows =
%       nodes)
%       element ==>  Point locations for each quadrilateral element
%
%     Date: August 2019 
%     Author: JT du Plessis

% -------------------------------------------------------------------------
% Init
% -------------------------------------------------------------------------

% Input arguments:
% rho = 0.5;
% run = 1;

%for vertices = 3:1:40
% num_vertices = 12; % Number of vertices
% Contour = [0 0 0; 1 0 0; 2 0 0; 3 0 0; 4 0 3];
%  Contour = [0 0 0; 1 0 0; 2 0 0; 3 0 0; 4 0 0.1; 5 0 0.2; 6 0.1 0.3; 7 0.2 0.4; 8 0.25 0.5; 9 0.2 0.6; 10 0.2 0.5; 11 0.2 0.4 ]; % Hardcode some contour
%   Contour = [0 0 0; 1 1 2; 2 2 4; 3 3 6];
%  Contour = [0 0 0; 0.5 0 0; 1 0 0; 1.5 0 0; 2 0 0];
%      Contour = [0 0 0; 1 1 1; 2 2 2; 3 3 3];
%  Contour = [0 0 0; 1 0 1; 2 0.5 1.5];
%   Contour = [0 0 0; 0.1 0 0; 0.2 0 0; 0.3 0 0; 0.4 0 0; 0.5 0 0;0.6 0 0; 0.7 0 0; 0.8 0 0; 0.9 0 0; 1 0 0; 1.1 0 0; 1.2 0 0; 1.3 0 0; 1.4 0 0; 1.5 0 0;1.6 0 0;1.7 0 0; 1.8 0 0;1.9 0 0; 2 0 0];
% Other:
num_nodes = length(Contour(:,1));
alpha = (2*pi) / num_vertices; % Angle step between vertices
points = zeros(num_vertices*(num_nodes-1),3);
seg = zeros(num_nodes-1,6);
displacement_vec = zeros(1,3);
norm_vec = zeros(1,3);

% -------------------------------------------------------------------------
% Main
% -------------------------------------------------------------------------

for node = 1:num_nodes
    % ----------------------------------
    % Form the points at the first node
    % ----------------------------------
    if node == 1
        
        norm_vec = Contour(node+1,:) - Contour(node,:); % Points away from next node
%         norm_vec(1:2) = norm_vec(1:2) +Contour(node,1:2);
        i = 2;
        basis_vec = [-norm_vec(mod(i,3)+1)/norm_vec(mod(i+1,3)+1) 1 0];
        while isnan(norm_vec(mod(i,3)+1)/norm_vec(mod(i+1,3)+1)) || isinf(norm_vec(mod(i,3)+1)/norm_vec(mod(i+1,3)+1)) || (i < 4)
            i = i + 1;
%             basis_vec = [ 0 1 -norm_vec(mod(i,3)+1)/norm_vec(mod(i+1,3)+1)];
            basis_vec = [-norm_vec(mod(i,3)+1)/norm_vec(mod(i+1,3)+1) 1 0];
%             basis_vec = [1 0 -norm_vec(mod(i,3)+1)/norm_vec(mod(i+1,3)+1)];
        end
%         if Contour(node,1) == 0
%            basis_vec = [ 0 -norm_vec(2)/norm_vec(1) 1];
% %         elseif Contour(node,3) == 0
%             basis_vec = [1 0 -norm_vec(2)/norm_vec(1)];
%         else
%             basis_vec = [-norm_vec(2)/norm_vec(1) 1 0]; % One possible basis vector orthogonal to norm_vec
%         end
%         basis_vec = basis_vec + Contour(node,:);
        basis_vec = rho* basis_vec/norm(basis_vec);
        norm_vec = norm_vec/norm(norm_vec); % Create unit vector
        
        for v = 1:num_vertices
            angle = alpha * (v-1);
            
            % Create set of points in plane by rotating point about normal
            % vector
            Q1 = [0 basis_vec];
            Q2 = [cos(angle/2),norm_vec(1)*sin(angle/2),norm_vec(2)*sin(angle/2),norm_vec(3)*sin(angle/2)];
            Q3 = quatmultiply(Q2,Q1);
            Q3 = quatmultiply(Q3,quatconj(Q2));
            points(v,1:3) = Q3(2:4);
            points(v,:) = points(v,:) + Contour(node,:); % Offset by starting node
%             if norm_vec(3) ~= 0
%                points(v,3) = Contour(node,3);
%             end

            
        end % Vertices
        
    else
        % ----------------------------------
        % If not the first node, then:
        %       1. Find plane that bisects the two normal vectors (vectors that 
        %          point to next and previous node)
        %       2. Project the points to the plane
        % ----------------------------------
        
        norm_vec = Contour(node-1,:) - Contour(node,:); % Points to previous node
        norm_vec = norm_vec/norm(norm_vec);
%         basis_vec = [-norm_vec(2)/norm_vec(1) 1 0];
        i = 2;
        if Contour(node,1) >= 5 && Contour(node,1) < 10
            basis_vec = [0 -norm_vec(mod(i,3)+1)/norm_vec(mod(i+1,3)+1) 1];
%         elseif Contour(node,1) > 10
%             basis_vec = [1 0 -norm_vec(mod(i,3)+1)/norm_vec(mod(i+1,3)+1)];
        else
            basis_vec = [0 -norm_vec(mod(i,3)+1)/norm_vec(mod(i+1,3)+1) 1];
        end
        while isnan(norm_vec(mod(i,3)+1)/norm_vec(mod(i+1,3)+1)) || isinf(norm_vec(mod(i,3)+1)/norm_vec(mod(i+1,3)+1))
            i = i + 1;
            basis_vec = [0 1 -norm_vec(mod(i,3)+1)/norm_vec(mod(i+1,3)+1)];
        end      
        
        if node == num_nodes
            basis_vec_rotated = norm_vec;
        else
            norm_vec_next = Contour(node+1,:) - Contour(node,:); % Points to next node
            % Calculate angle between vector pointing to previous node and
            % vector pointing to next node
            angle_between = abs(acos(dot(norm_vec_next,norm_vec)/(norm(norm_vec_next)*norm(norm_vec))));
            if angle_between == pi || angle_between == 0
                bisector = basis_vec;
            else
                bisector = norm_vec/norm(norm_vec) + norm_vec_next/norm(norm_vec_next);
%                 bisector = norm_vec + norm_vec_next;
            end

            bisector = bisector/norm(bisector);
            basis_vec2 = cross(bisector,norm_vec);
            basis_vec2 = basis_vec2/norm(basis_vec2);
                 
            % Rotate bisector 90 degrees about basis_vec2 to form normal
            % vector to tilted plane
            beta = pi/2;
            Q1 = [0 bisector];
            Q2 = [cos(beta/2),basis_vec2(1)*sin(beta/2),basis_vec2(2)*sin(beta/2),basis_vec2(3)*sin(beta/2)];
            Q3 = quatmultiply(Q2,Q1);
            Q3 = quatmultiply(Q3,quatconj(Q2));
            basis_vec_rotated = Q3(2:4); % This vector is normal to tilted plane
            % basis_vec_rotated should equal norm_vec when beta = pi
            
        end

        % Project points to plane
        for v = 1:num_vertices
             
            point_index      = v + num_vertices*(node-1);
            point_index_prev = v + num_vertices*(node-2);
            
            vec1 = points(point_index_prev,:) - Contour(node,:); % 
            dist = dot(vec1,basis_vec_rotated); % Perpendicular distance to plane
            
            if dist < 0
                dist = -1 * dist;
            end
            % Translate point to plane, in direction of norm_vec
            points(point_index,1:3) =  points(point_index_prev,:) - (dist*norm_vec);
        end

        
    end
end

% This ensures the first and last nodes are on the plate
% points(1:num_vertices,3) = Contour(1,3);
% points(end-num_vertices+1:end,3) = Contour(end,3);

% ------------------------------------------------------------------------
% -------------------------New 13/05/2020---------------------------------
% ------------------------------------------------------------------------
% Mesh connection to external structure
% Add points for connection to external plate/structure

% if cyl_def.firstNode == "conn" || cyl_def.lastNode == "conn" 
%    point_index      = point_index - num_vertices;               % Start at beginning of last node
% %    point_index_next = mod(point_index+1, num_vertices) + p1 - num_vertices;    % Next point/vertex
%    point_index_next = mod(point_index, num_vertices) + 1;    % Next point/vertex
%     for p = 1:num_vertices
%         connection_index = point_index + p + num_vertices; % Create new set of points
%         % Create points at first node
%         points(connection_index,1:3) = 2 * cross((points(point_index+p,:) + points(point_index_next+p,:))/2,norm_vec); % Scaled point in midpoint of prev and next point
%     end
%     % Repeat elements for last node
%     points(end+1:end + num_vertices,:) = [repelem(Contour(end,1),num_vertices)',points(end - (num_vertices) +1:end,2:3)]; 
%     
% end




% Now take the points and form the quadrilateral elements
element = zeros((num_nodes-2)*(num_vertices)-1,3);
for p = 1:(num_nodes-1)*(num_vertices)
    for corner = 1:4
        %1,2,7,8
        %2,3,8,9
        offset_curr_node = abs(1-mod(corner,2));
        offset_next_node = 0;
        if corner > 2
            offset_next_node = num_vertices;
        end     
        if (corner == 2 || corner == 4) && mod((p+offset_curr_node),num_vertices) == 1
            offset_curr_node = -(num_vertices-1); %go back to first vertice in node
        end      
        element(p,corner) = p+offset_curr_node+offset_next_node;
    end
end

% Swap last nodes to get sequential indices going counter-clockwise
third = element(:,3);
fourth = element(:,4);
element(:,3) = fourth;
element(:,4) = third;

offset = 0;

% ------------------------------------------------------------------------
% -------------------------New 13/05/2020---------------------------------
% ------------------------------------------------------------------------
% if cyl_def.firstNode == "conn"
%     offset = num_vertices*4;
%     for v = 1:num_vertices % There will be num_vertices extra points at first and last node
%         % The elements at the start
%         % 'p' is the last element (#elements)
%         element(p+v, 1) = element(v,1);
%         element(p+v, 2) = element(v,2);
%         element(p+v, 3) = element(p,4) + v; % Add extra
%         element(p+v, 4) = element(p,4) + v; % Degenerate coordinate
%         
%         % The elements at the end
%         element(p+v+num_vertices, 1) = element(p-num_vertices+v,4);
%         element(p+v+num_vertices, 2) = element(p-num_vertices+v,3);
%         element(p+v+num_vertices, 3) = element(p,4) + v + num_vertices; % Add extra
%         element(p+v+num_vertices, 4) = element(p,4) + v + num_vertices; % Degenerate coordinate
%     end
%     p = length(element);
%     element(p-(2*num_vertices)+1:p-num_vertices,3:4) = circshift(element(p-(2*num_vertices)+1:p-num_vertices,3:4),round(num_vertices/4));
%     element(p-(num_vertices)+1:p,3:4) = circshift(element(p-(num_vertices)+1:p,3:4),round(num_vertices/4));
%     
%     
% %     % Now to form secondary elements between each 'flag' element
% %     for v = 1:num_vertices
% %         next_ind = p- (2*num_vertices) + mod(v-1,num_vertices)+1;
% %         prev_ind = p- (2*num_vertices) + mod(v-2,num_vertices)+1;
% %        % Elements at the start node
% %        element(p+v,1) = element(v,1);
% %        element(p+v,2) = element(next_ind,3);
% %        element(p+v,3) = element(prev_ind,3);
% %        element(p+v,4) = element(prev_ind,3); % Degenerate. Don't really care about this.
% %        
% %        % Elements at the last node
% %        element(p+v+num_vertices,1) = element(p-num_vertices+v,1);
% %        element(p+v+num_vertices,2) = element(next_ind + num_vertices,3);
% %        element(p+v+num_vertices,3) = element(prev_ind + num_vertices,3);
% %        element(p+v+num_vertices,4) = element(prev_ind + num_vertices,3); % Degenerate. Don't really care about this.
%        
%        
% %     end
% %    
% %     TR = triangulation(tri_nodes,points);
% %     N = nearestNeighbor(TR, points);
% %     testd = delaunayTriangulation(points);
% end
%  [tri_nodes,triangles] = QuadtoTri(element, num_vertices, cyl_def);
% 
%  PlotMesh(points,tri_nodes, tri_nodes, 1);
% ------------------------------------------------------------------------
% -------------------------New 19/02/2020---------------------------------
% ------------------------------------------------------------------------
% Mesh degenerate quadrilaterals on the endcaps
if cyl_def.firstNode == "endCap" || cyl_def.lastNode == "endCap"
    p = length(element);
    point_index = length(points);
    
    if cyl_def.firstNode == "endCap" && cyl_def.lastNode == "endCap"
        first_index = point_index +1;
        second_index = point_index +2;
        second_offset = num_vertices;
        points(first_index, 1:3) = [0,0,0];
        points(second_index, 1:3) = [0,0,0];
    elseif cyl_def.lastNode == "endCap" && cyl_def.firstNode ~= "endCap"
       second_offset = 0;
        first_index = 0;
        second_index = point_index +1;
        points(second_index, 1:3) = [0,0,0];
    elseif cyl_def.lastNode ~= "endCap" && cyl_def.firstNode == "endCap"
        second_offset = 0;
        first_index = point_index +1;
        second_index = 0;
        points(first_index, 1:3) = [0,0,0];
    end
    
    for v = 1:num_vertices % There will be num_vertices degenerate quads on the endcaps
        
        if cyl_def.firstNode == "endCap"
            % Average out the points at the end nodes
             points(first_index, 1:3) = points(first_index, 1:3) + points(v,1:3); % First node
            % The elements at the start
            % 'p' is the last element (#elements)
            element(p+v, 1) = element(v,1);
            element(p+v, 2) = element(v,2);
            element(p+v, 3) = first_index; % Second to last point in points matrix
            element(p+v, 4) = first_index;
        end

        % The elements at the end
        if cyl_def.lastNode == "endCap"
             % Average out the points at the end nodes
            points(second_index, 1:3) = points(second_index, 1:3) + points(point_index-v+1,1:3); % Last node
            
            element(p+v+second_offset, 1) = element(p-offset-v+1,3);
            element(p+v+second_offset, 2) = element(p-offset-v+1,4);
            element(p+v+second_offset, 3) = second_index; % Last point in points matrix
            element(p+v+second_offset, 4) = second_index;
        end
    end
    
    if cyl_def.firstNode == "endCap"
         points(first_index,1:3) = points(first_index,1:3)/num_vertices;
    end
    if cyl_def.lastNode == "endCap"
         points(second_index,1:3) = points(second_index,1:3)/num_vertices;
    end
   
end




quadElements = {};
for node = 1:length(element)
    quadElements = [quadElements;element(node,1:4),[0 0 0 0 0 0 0 0]]; 
end


% [indexes,error(run)] = TestPlanar(points,element);
% [new_triangles,triangles] = QuadtoTri(element);
% PlotMesh(points,element,triangles);

% run = run + 1;
% end
% vertices = 3:1:40;
% figure();
% plot(vertices,error);
% title('Percentage Error vs Number of Vertices');
% xlabel('Vertices');
% ylabel('Percentage Error (%)');
 end
