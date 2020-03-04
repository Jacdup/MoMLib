 function [points,quadElements, element] = QuadMesh_v5(Contour,num_vertices,rho)
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
run = 1;

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
        
        basis_vec = [-norm_vec(2)/norm_vec(1) 1 0]; % One possible basis vector orthogonal to norm_vec
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
            points(v,2:3) = Q3(3:4);
            points(v,1) = Contour(node,1);
            
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
        basis_vec = [-norm_vec(2)/norm_vec(1) 1 0];
        
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
            point_index = v + num_vertices*(node-1);
            point_index_prev = v + num_vertices*(node-2);
            
            vec1 = points(point_index_prev,:) - Contour(node,:);
            dist = dot(vec1,basis_vec_rotated);
            
            if dist < 0
                dist = -1 * dist;
            end
            % Translate point to plane, in direction of norm_vec
            points(point_index,1:3) =  points(point_index_prev,:) - (dist*norm_vec);
        end
        
    end
end


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

% ------------------------------------------------------------------------
% -------------------------New 19/02/2020---------------------------------
% ------------------------------------------------------------------------
% Mesh degenerate quadrilaterals on the endcaps
% points(point_index + 1, 1:3) = Contour(1,:);
% points(point_index + 2, 1:3) = Contour(end,:);
% for v = 1:num_vertices % There will be num_vertices degenerate quads on the endcaps
%     % The elements at the start
%     % 'p' is the last element (#elements)
%     element(p+v, 1) = element(v,1);
%     element(p+v, 2) = element(v,2);
%     element(p+v, 3) = element(p,4) + 1; % Add extra
%     element(p+v, 4) = element(p,4) + 1; % Degenerate coordinate
%     
%     % The elements at the end
%     element(p+v+num_vertices, 1) = element(p-v+1,4);
%     element(p+v+num_vertices, 2) = element(p-v+1,3);
%     element(p+v+num_vertices, 3) = element(p,4) + 2; % Add extra
%     element(p+v+num_vertices, 4) = element(p,4) + 2; % Degenerate coordinate
% end


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
