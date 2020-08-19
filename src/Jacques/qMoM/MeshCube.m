function [node_coords, quad_elems] = MeshCube(points, radius)

a = radius; % Radius
increment = 2/(points-1);
x = -1:increment:1;
y = -1:increment:1;
z = -1:increment:1;
[X,Y,Z] = meshgrid(x);
points = [X(:),Y(:),Z(:);zeros(1,3)];
new_index = 1;
points_per_dim = length(z(1,:,:));


% r_s = zeros(length(points(:,1)),3);
for i = 1:length(points(:,1))
    % Remove non-face entries
    if any(points(i,:) == -1) || any(points(i,:) == 1)
        new_Points(new_index, :) = points(i,:);
        new_index = new_index + 1;
    end
end
% Sort according to 
% B = sortrows(new_Points, [1 2 3]);
% scatter3(new_Points(:,1),new_Points(:,2),new_Points(:,3))
% patch(new_Points(:,1),new_Points(:,2),new_Points(:,3));

% Form quadrilateral elements
% for p = 1:8 % All faces
points_per_face = points_per_dim*(points_per_dim);
iter = 0;
iter1 = 0;
left_offset = 0;
    for f = 1:points_per_face  % All points in face       
            if mod(f+iter,points_per_dim) == 0
                iter = iter + 1;
                iter1 = iter1 + 1;
            end
    end


elem_idx = 0;
% face1_idx = 1:points_per_face;
% face2_idx = points_per_face:(points_per_face*2);
% face3_idx = points_per_face:(points_per_face*3);
% face4_idx = points_per_face:(points_per_face*4);
% face5_idx = points_per_face:(points_per_face*5);
quad_element = zeros(6,points_per_face);
face1_idx = 1;
face2_idx = points_per_face;
face3_idx = face2_idx + points_per_face;
node = 1;
iter = 0;
iter2 = 0;
iter3 = 0;
iter4 = 0;
iter5 = 0;
iter6 = 0;
for ii = 1:length(new_Points(:,1))
    
    % Bottom face
    if (new_Points(ii,3) == -1)
        iter = iter + 1;
        quad_element(1,iter) = ii; % Need to reshape this later
    end
    
    % Left face
    if (new_Points(ii,1) == -1)
        iter2 = iter2 + 1;
        quad_element(2,iter2) = ii;
    end
        % Front face
    if (new_Points(ii,2) == -1)
        iter3 = iter3 + 1;
        quad_element(3,iter3) = ii;
    end
      if (new_Points(ii,1) == 1)
        iter4 = iter4 + 1;
        quad_element(4,iter4) = ii;
      end
      if (new_Points(ii,2) == 1)
        iter5 = iter5 + 1;
        quad_element(5,iter5) = ii;
      end
      if (new_Points(ii,3) == 1)
        iter6 = iter6 + 1;
        quad_element(6,iter6) = ii;
    end
    
end

for i = 1:points_per_dim
     % Selects every 6 elements
    elems_f1(i,:) = quad_element(1,i:points_per_dim:end);
    elems_f2(i,:) = quad_element(2,i:points_per_dim:end);
    elems_f3(i,:) = quad_element(3,i:points_per_dim:end); 
    elems_f4(i,:) = quad_element(4,i:points_per_dim:end);
    elems_f5(i,:) = quad_element(5,i:points_per_dim:end); 
    elems_f6(i,:) = quad_element(6,i:points_per_dim:end);
end
iter = 1;
% elems_f1_new = zeros(25,4);
% elems_f2_new = zeros(25,4);
% elems_f3_new = zeros(25,4);
% elems_f4_new = zeros(25,4);
% elems_f5_new = zeros(25,4);
% elems_f6_new = zeros(25,4);
%Select every 'block' in the matrix to form the elements
for i = 1:points_per_dim-1
    for j = 1:points_per_dim-1
         elems_f1_new(iter,:) = reshape(elems_f1(i:i+1,j:j+1),[1,4]); 
         elems_f2_new(iter,:) = reshape(elems_f2(i:i+1,j:j+1),[1,4]); 
         elems_f3_new(iter,:) = reshape(elems_f3(i:i+1,j:j+1),[1,4]); 
         elems_f4_new(iter,:) = reshape(elems_f4(i:i+1,j:j+1),[1,4]); 
         elems_f5_new(iter,:) = reshape(elems_f5(i:i+1,j:j+1),[1,4]); 
         elems_f6_new(iter,:) = reshape(elems_f6(i:i+1,j:j+1),[1,4]); 
         iter = iter + 1;
    end
end
% elems_f1_reshape = reshape(elems_f1(:,:),[9, 4]);
% elems = reshape(quad_element(:,:),[108, 2]);
elems = [elems_f1_new;elems_f2_new;elems_f3_new;elems_f4_new;elems_f5_new;elems_f6_new];

% Swap last nodes to get sequential indices going counter-clockwise
elems_f1_new(:,[3 4]) = elems_f1_new(:,[4 3]);
elems_f2_new(:,[3 4]) = elems_f2_new(:,[4 3]);
elems_f3_new(:,[3 4]) = elems_f3_new(:,[4 3]);
elems_f4_new(:,[3 4]) = elems_f4_new(:,[4 3]);
elems_f5_new(:,[3 4]) = elems_f5_new(:,[4 3]);
elems_f6_new(:,[3 4]) = elems_f6_new(:,[4 3]);
elems(:,[3 4]) = elems(:,[4 3]);

% Make sphere
for new_index = 1:length(new_Points(:,1))
        x = new_Points(new_index,1);
        y = new_Points(new_index,2);
        z = new_Points(new_index,3);
        
        r_s(new_index,1) = x*sqrt(1-((y^2)/2) - ((z^2)/2) + (((y^2)*(z^2))/3));
        r_s(new_index,2) = y*sqrt(1-((z^2)/2) - ((x^2)/2) + (((z^2)*(x^2))/3));
        r_s(new_index,3) = z*sqrt(1-((x^2)/2) - ((y^2)/2) + (((x^2)*(y^2))/3));

%         r_s(new_index,1) = new_Points(new_index,1)/norm(new_Points(new_index,1:3));
%         r_s(new_index,2) = new_Points(new_index,2)/norm(new_Points(new_index,1:3));
%         r_s(new_index,3) = new_Points(new_index,3)/norm(new_Points(new_index,1:3));
        r_s(new_index,:) = r_s(new_index,:) * a;
end
% patch('Faces', elems(:,:), 'Vertices',new_Points, 'FaceColor','Green' )
c1 = repelem(0.1,points_per_dim)';
c2 = repelem(0.5,points_per_dim)';
col = [c1; c2;c1;c2;c1;c2];
% figure
% patch('Faces', elems_f1_new(:,:), 'Vertices',r_s, 'FaceColor', 'r')
% patch('Faces', elems_f2_new(:,:), 'Vertices',r_s, 'FaceColor', 'b')
% patch('Faces', elems_f3_new(:,:), 'Vertices',r_s, 'FaceColor', 'g')
% patch('Faces', elems_f4_new(:,:), 'Vertices',r_s, 'FaceColor', 'y')
% patch('Faces', elems_f5_new(:,:), 'Vertices',r_s, 'FaceColor', 'm')
% patch('Faces', elems_f6_new(:,:), 'Vertices',r_s, 'FaceColor', 'c')
% axis equal

% figure
% patch('Faces', elems_f1_new(:,:), 'Vertices',new_Points, 'FaceColor', 'r')
% patch('Faces', elems_f2_new(:,:), 'Vertices',new_Points, 'FaceColor', 'b')
% patch('Faces', elems_f3_new(:,:), 'Vertices',new_Points, 'FaceColor', 'g')
% patch('Faces', elems_f4_new(:,:), 'Vertices',new_Points, 'FaceColor', 'y')
% patch('Faces', elems_f5_new(:,:), 'Vertices',new_Points, 'FaceColor', 'm')
% patch('Faces', elems_f6_new(:,:), 'Vertices',new_Points, 'FaceColor', 'c')
% axis equal

node_coords = r_s;
quad_elems = elems;
% patch('Faces', elems(:,:), 'Vertices',r_s, 'FaceColor','Green' )
% plot3(r_s(:,1),r_s(:,2),r_s(:,3))
% patch(r_s(:,1),r_s(:,2),r_s(:,3));
% view(3);