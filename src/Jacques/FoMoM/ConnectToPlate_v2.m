function [node_coords, elements, plate_tri_nodes, cyl_def, last_element_val, plate_node_coords] = ConnectToPlate_v2(cyl_def, node_coords,elements, vertices)
% All this function should do, is search the FEKO Nastran mesh for the
% polygon_points (connection points), and replace the cylinder end nodes by 
% those FEKO nodes which correspond to the coordinates in the FEKO mesh.


% if isempty(polygon_points)
    polygon_points_1 = node_coords(1:vertices,:);
    polygon_points_2 = node_coords(end-vertices+1:end,:);
    
% else
%     polygon_points_1 = polygon_points(1:vertices,:);
%     polygon_points_2 = polygon_points(vertices+1:end,:);
% end

% polygon_points_3 = polygon_points_1; % VERY hardcoded here for y-axis offset of -1
% polygon_points_3(:,2) = polygon_points_3(:,2) - 1;% Shift y-axis by 1
% polygon_points_4 = polygon_points_2; % VERY hardcoded here for y-axis offset of -1
% polygon_points_4(:,2) = polygon_points_4(:,2) - 1;% Shift y-axis by 1


%  polygon_points_1 = [0.500000000000000	0.650000000000000	0
%         0.453647450843758	0.642658477444273	0
%         0.411832212156129	0.621352549156242	0
%         0.378647450843758	0.588167787843871	0
%         0.357341522555727	0.546352549156242	0
%         0.350000000000000	0.500000000000000	0
%         0.357341522555727	0.453647450843758	0
%         0.378647450843758	0.411832212156129	0
%         0.411832212156129	0.378647450843758	0
%         0.453647450843758	0.357341522555727	0
%         0.500000000000000	0.350000000000000	0
%         0.546352549156242	0.357341522555727	0
%         0.588167787843871	0.378647450843758	0
%         0.621352549156242	0.411832212156129	0
%         0.642658477444273	0.453647450843758	0
%         0.650000000000000	0.500000000000000	0
%         0.642658477444273	0.546352549156242	0
%         0.621352549156242	0.588167787843871	0
%         0.588167787843871	0.621352549156242	0
%         0.546352549156242	0.642658477444273	0];

%      polygon_points_2 = polygon_points_1;
%     polygon_points_2(:,3) = polygon_points_2(:,3) + 3.0;



%          [plate_node_coords, plate_tri_nodes] = ImportTriangleMeshNastran('C:\Users\19083688\Desktop\Masters\FEKO Models\Meshes\square_plate_hole_regular.nas');
if cyl_def.firstNode == "conn" && cyl_def.lastNode ~= "conn"
%                 [plate_node_coords, plate_tri_nodes] = ImportTriangleMeshNastran('C:\Users\19083688\Desktop\Masters\FEKO Models\Meshes\2_square_plate_hole_20v_fine.nas');
    [plate_node_coords, plate_tri_nodes] = ImportTriangleMeshNastran('C:\Users\19083688\Desktop\Masters\FEKO Models\Meshes\square_plate_hole_20v_regular.nas');
else
%      [plate_node_coords, plate_tri_nodes] = ImportTriangleMeshNastran('C:\Users\19083688\Desktop\Masters\FEKO Models\Meshes\2_square_plate_hole_20v_fine.nas');
%                [plate_node_coords, plate_tri_nodes] = ImportTriangleMeshNastran('C:\Users\19083688\Desktop\Masters\FEKO Models\Meshes\square_plate_hole_20v_regular.nas');
    %            [plate_node_coords, plate_tri_nodes] = ImportTriangleMeshNastran('C:\Users\19083688\Desktop\Masters\FEKO Models\Meshes\square_plate_2holes.nas');
    %             [plate_node_coords, plate_tri_nodes] = ImportTriangleMeshNastran('C:\Users\19083688\Desktop\Masters\FEKO Models\Meshes\square_plate_2holes_smaller.nas');
    [plate_node_coords, plate_tri_nodes] = ImportTriangleMeshNastran('C:\Users\19083688\Desktop\Masters\FEKO Models\Meshes\square_plate_4holes.nas');
    
end



if cyl_def.firstNode == "conn" && cyl_def.lastNode ~= "conn"
    [~, LocB] = ismembertol(polygon_points_1, plate_node_coords, 0.001, 'ByRows', 1); % Determine where the polygon points are in the FEKO mesh
    plate_polygon_nodes_1 = nonzeros(LocB);
    
elseif cyl_def.firstNode == "conn" && cyl_def.lastNode == "conn"
    [~, LocB] = ismembertol(polygon_points_1, plate_node_coords, 0.001, 'ByRows', 1); % Determine where the polygon points are in the FEKO mesh
    plate_polygon_nodes_1 = nonzeros(LocB);
    [~, LocB] = ismembertol(polygon_points_2, plate_node_coords, 0.001, 'ByRows', 1); % Determine where the polygon points are in the FEKO mesh
    plate_polygon_nodes_2 = nonzeros(LocB);
    
    
%     if cyl_def.coupling % Then there are 4 plate connections
%         [~, LocB] = ismembertol(polygon_points_3, plate_node_coords, 0.001, 'ByRows', 1); % Determine where the polygon points are in the FEKO mesh
%         plate_polygon_nodes_3 = nonzeros(LocB);
%         [~, LocB] = ismembertol(polygon_points_4, plate_node_coords, 0.001, 'ByRows', 1); % Determine where the polygon points are in the FEKO mesh
%         plate_polygon_nodes_4 = nonzeros(LocB);
%     end
    
    %             elements = elements - (2*vertices);
    %             last_element_val = max(max(elements));
else
    [~, LocB] = ismembertol(polygon_points_2, plate_node_coords, 0.001, 'ByRows', 1); % Determine where the polygon points are in the FEKO mesh
    plate_polygon_nodes_2 = nonzeros(LocB);
    
end


% ---------------------------------------------------------------------------------------------
% Method Number 1
% ---------------------------------------------------------------------------------------------
%         [node_coords1, polygon_cyl_elements] = ConnectToPlate(polygon_points_1,plate_polygon_nodes_1, node_coords, elements, vertices, cyl_def);
%          new_elements = [plate_polygon_nodes,];
% ---------------------------------------------------------------------------------------------
% Method Number 2
% ---------------------------------------------------------------------------------------------
%         %           plate_node_coords(LocB,:) = node_coords(1:vertices,:); % These two are the same
%         %           plate_node_coords(LocB,:) = []; % Remove all the FEKO polygon points
%         %           plate_tri_nodes = plate_tri_nodes + vertices;
%         last_element_val = max(max(elements));
%         temp = 1:vertices;
%         LocB = LocB + last_element_val;
%         plate_tri_nodes = plate_tri_nodes + last_element_val;
%
% %         plate_node_coords(LocB,:) = [];
%
%         for i = 1:3 % Replace all the FEKO polygon nodes by our mesh nodes
%             Out = {};
%             for iA = 1:numel(LocB)
%                 Out{iA} = find(plate_tri_nodes(:,i) == LocB(iA));
%                 if Out{iA} ~= 0
%                     plate_tri_nodes(Out{iA},i) = temp(iA);
%                 end
%             end
%         end
%
%         clear temp LocB1 LocB2;
%         node_coords = [node_coords; plate_node_coords];
%         %          elements = [elements;polygon_cyl_elements];
%
%         cyl_def.num_plate_nodes = length(plate_tri_nodes) + 1; % Number of nodes to exclude from MBF
%         cyl_def.plate_polygon_nodes = plate_polygon_nodes;
%
%
%         [tri_nodes,triangles] = QuadtoTri(elements, vertices,cyl_def);
%         tri_nodes = [tri_nodes;(plate_tri_nodes)];
% ---------------------------------------------------------------------------------------------
% Method Number 3
% ---------------------------------------------------------------------------------------------

if cyl_def.firstNode == "conn" && cyl_def.lastNode ~= "conn" % First node is connection
    elements = elements - vertices;
    last_element_val = max(max(elements));
    
    elements(1:vertices,1) = plate_polygon_nodes_1 + last_element_val;
    elements(1:vertices,2) = circshift(plate_polygon_nodes_1,-1) + last_element_val;
    node_coords(1:vertices,:) = [];
    
    cyl_def.plate_polygon_nodes = plate_polygon_nodes_1;
    
elseif cyl_def.firstNode == "conn" && cyl_def.lastNode == "conn"            % Both are connections
    elements = elements - vertices;                                         % Offset all cylinder element indices (since they will be overwritten by the FEKO nodes)
    elements(end-vertices+1:end,:) = elements(end-vertices+1:end,:) - vertices; % Offset last (vertices) element indices
    last_element_val = max(max(elements));                                  % Highest triangle number in cylinder mesh
    
    elements(1:vertices,1) = plate_polygon_nodes_1 + last_element_val;      % Add FEKO element numbers of first node (at first connection), offset by highest triangle number in cylinder mesh
    elements(1:vertices,2) = circshift(plate_polygon_nodes_1,-1) + last_element_val; %Add FEKO element numbers of second triangle node
    node_coords(1:vertices,:) = [];                                         % Remove redundant coordinates of cylinder
    
    elements(end-vertices+1:end,1) = plate_polygon_nodes_2 +  last_element_val; % Add FEKO element numbers of first node, second connection
    elements(end-vertices+1:end,2) = circshift(plate_polygon_nodes_2,-1) +  last_element_val; %Add FEKO element numbers of second node, second connection
    node_coords(end-vertices+1:end,:) = [];                                 % Remove redundant coordinates of cylinder
    
    cyl_def.plate_polygon_nodes = [cyl_def.plate_polygon_nodes;plate_polygon_nodes_1];                    % Add the nodes to the struct for future reference
    cyl_def.plate_polygon_nodes_end = [cyl_def.plate_polygon_nodes_end;plate_polygon_nodes_2];
    cyl_def.last_element_val = last_element_val;
    
    
else % Last node is a connection
    elements(end-vertices+1:end,:) = elements(end-vertices+1:end,:) - vertices;
    last_element_val = max(max(elements));
    
    elements(end-vertices+1:end,1) = plate_polygon_nodes_2 + last_element_val;
    elements(end-vertices+1:end,2) = circshift(plate_polygon_nodes_2,-1) + last_element_val;
    node_coords(end-vertices+1:end,:) = [];
    
    cyl_def.plate_polygon_nodes_end = plate_polygon_nodes_2;
end

if cyl_def.coupling == 0
     node_coords = [node_coords; plate_node_coords];
end
%          elements = [elements;polygon_cyl_elements];

cyl_def.num_plate_nodes = length(plate_tri_nodes) + 1; % Number of nodes to exclude from MBF




end
