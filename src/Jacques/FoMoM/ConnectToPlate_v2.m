function [node_coords, elements, plate_tri_nodes, cyl_def, last_element_val] = ConnectToPlate_v2(cyl_def, node_coords,elements, vertices)

polygon_points_1 = node_coords(1:vertices,:);
polygon_points_2 = node_coords(end-vertices+1:end,:);
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
%             [plate_node_coords, plate_tri_nodes] = ImportTriangleMeshNastran('C:\Users\19083688\Desktop\Masters\FEKO Models\Meshes\2_square_plate_hole_20v_fine.nas');
       [plate_node_coords, plate_tri_nodes] = ImportTriangleMeshNastran('C:\Users\19083688\Desktop\Masters\FEKO Models\Meshes\square_plate_hole_20v_regular.nas');
       else
%            [plate_node_coords, plate_tri_nodes] = ImportTriangleMeshNastran('C:\Users\19083688\Desktop\Masters\FEKO Models\Meshes\square_plate_hole_20v_regular.nas');
%            [plate_node_coords, plate_tri_nodes] = ImportTriangleMeshNastran('C:\Users\19083688\Desktop\Masters\FEKO Models\Meshes\square_plate_2holes.nas');
            [plate_node_coords, plate_tri_nodes] = ImportTriangleMeshNastran('C:\Users\19083688\Desktop\Masters\FEKO Models\Meshes\square_plate_2holes_smaller.nas');
       end

        if cyl_def.firstNode == "conn" && cyl_def.lastNode ~= "conn"
            [~, LocB] = ismembertol(polygon_points_1, plate_node_coords, 0.001, 'ByRows', 1); % Determine where the polygon points are in the FEKO mesh
            plate_polygon_nodes_1 = nonzeros(LocB);
            
        elseif cyl_def.firstNode == "conn" && cyl_def.lastNode == "conn"
            [~, LocB] = ismembertol(polygon_points_1, plate_node_coords, 0.001, 'ByRows', 1); % Determine where the polygon points are in the FEKO mesh
            plate_polygon_nodes_1 = nonzeros(LocB);
            [~, LocB] = ismembertol(polygon_points_2, plate_node_coords, 0.001, 'ByRows', 1); % Determine where the polygon points are in the FEKO mesh
            plate_polygon_nodes_2 = nonzeros(LocB);
            
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
        cyl_def.plate_polygon_nodes_end = [];

        if cyl_def.firstNode == "conn" && cyl_def.lastNode ~= "conn" % First node is connection
            elements = elements - vertices;
            last_element_val = max(max(elements));
            
            elements(1:vertices,1) = plate_polygon_nodes_1 + last_element_val;
            elements(1:vertices,2) = circshift(plate_polygon_nodes_1,-1) + last_element_val;
            node_coords(1:vertices,:) = [];
            
            cyl_def.plate_polygon_nodes = plate_polygon_nodes_1;
            
        elseif cyl_def.firstNode == "conn" && cyl_def.lastNode == "conn" % Both are connections
            elements = elements - vertices;
            elements(end-vertices+1:end,:) = elements(end-vertices+1:end,:) - vertices;
            last_element_val = max(max(elements));
            cyl_def.last_element_val = last_element_val;
%             if length(plate_polygon_nodes_1) < vertices
%                 ""
%             end
            
            elements(1:vertices,1) = plate_polygon_nodes_1 + last_element_val;
            elements(1:vertices,2) = circshift(plate_polygon_nodes_1,-1) + last_element_val;
            node_coords(1:vertices,:) = [];
            
            elements(end-vertices+1:end,1) = plate_polygon_nodes_2 +  last_element_val;
            elements(end-vertices+1:end,2) = circshift(plate_polygon_nodes_2,-1) +  last_element_val;
            node_coords(end-vertices+1:end,:) = [];
            
            cyl_def.plate_polygon_nodes = plate_polygon_nodes_1;
            cyl_def.plate_polygon_nodes_end = plate_polygon_nodes_2;
            
        else % Last node is a connection
            elements(end-vertices+1:end,:) = elements(end-vertices+1:end,:) - vertices;
            last_element_val = max(max(elements));
            
            elements(end-vertices+1:end,1) = plate_polygon_nodes_2 + last_element_val;
            elements(end-vertices+1:end,2) = circshift(plate_polygon_nodes_2,-1) + last_element_val;
            node_coords(end-vertices+1:end,:) = [];
            
            cyl_def.plate_polygon_nodes_end = plate_polygon_nodes_2;
        end  
        
        node_coords = [node_coords; plate_node_coords];
        %          elements = [elements;polygon_cyl_elements];
        
        cyl_def.num_plate_nodes = length(plate_tri_nodes) + 1; % Number of nodes to exclude from MBF
        

        
        
end
