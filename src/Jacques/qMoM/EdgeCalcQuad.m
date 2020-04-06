function [quadElements,N, basis_supports] = EdgeCalcQuad(quadElements,vertices)%filename recieved is a string
%Jacques T du Plessis
%September 2019
%19083688@sun.ac.za
%
% Extension to Robey's code for quadrilateral support

%DESCRIPTION
%This routine allocates the current directions over the entire mesh. It
%must be done flawlessly every time. The current is returned using the
%current direction is defined over the edge as A B C D a b c d, where 'a'
%determines the current direction over the BC '1' means the current flows
%out and '-1' means the current flows in.
%
%This functions allocates the current directions and DOFF numbering over
%the quadrilateral mesh. It recieves a sorted quad cell matrix.
%The layout is :
%       For quadrilaterals:
%       A  B  C  D || E  F  G  H || I  J  K  M
%       
%       current directions:
%       (+1) current going out (-1) coming in
%       E -> current direction over edge AC
%       F -> current direction over edge CD
%       G -> current direction over edge BD
%       H -> current direction over edge AB
%
%       DOF numbering:
%       if (-1) not a DOF otherwise DOF number
%       I (9) -> number for edge AC
%       J (10)-> number for edge CD
%       K (11)-> number for edge BD
%       M (12)-> number for edge AB
%
%UPDATE 24/04/2019
%The code is updated from version 3.0 to now also return the basis function
%supports, an N by 2 array that contains the two triangles supporting the
%basis function at the nth DOFF(the row index is the doff number it supports)
%This is a trivial implementation.

%Define variables and constants
connectivity_data = {};
% edge_select = [1 2 ;1 3 ;2 3]; AB; AC; BC
 
%  edge_select = [1 2; 2 3; 4 3; 1 4];
edge_select = [1 2;1 4 ; 4 3; 2 3];
% Edge_select = [A B; B C; D C; A D];
% direction_select = [6 5 4];
direction_select = [8 5 6 7];
% doff_select = [9 8 7];
doff_select = [12 9 10 11];
DOFF_NUM = 0;
% tri = [0, 0, 0];
% quad = [0, 0, 0, 0];
high_num = 0;
N = 0;
basis_supports = [];

%populate connectivity list data
for quad_index = 1:size(quadElements,1)
    quad = quadElements{quad_index,1:4};
    for node = 1:4
      
        if high_num >= quad(node)%means we dont have to add another cell to the end
            
            connectivity_data{quad(node)} = [connectivity_data{quad(node)};quad_index];
        else %new node, create new cell
            
            connectivity_data{quad(node)} = [quad_index];
            high_num = quad(node);
        end
    end
end
%  connectivity_data{1}
% connectivity_data{2}
% connectivity_data{3}
% connectivity_data{4}
% connectivity_data{5}
% connectivity_data{6}


%-----------%-----------%-----------%-----------%-----------%-----------%-----------%-----------
%-----------Now we will simultaneously number the DOF's and assign the current-----------
%-----------directions for the quadrilaterals.-----------
%-----------%-----------%-----------%-----------%-----------%-----------%-----------%-----------


%Now we will iterate through the quadrilaterals

for quad_index = 1:size(quadElements,1)
    quad = quadElements{quad_index,1:4};
    for counter = 1:4
        %once we have selected a quadrilateral T_n then we will select an edge
        %from that quadrilateral
        curr_edge = [quad(edge_select(counter,1)),quad(edge_select(counter,2))];
        %once the edge of the quadrilateral is selected then we go to the
        %connectivity_list and check if the edge shares two quadElements, if
        %the edges does indeed share two elements then it is a DOF. see
        %below.
        
        %choose the relevant connectivity list columns
        node_A = connectivity_data{curr_edge(1)};
        node_B = connectivity_data{curr_edge(2)};
        common_quads = intersect(node_A,node_B);
        %size(common_triangles,1)
        %If the length of common_triangles is equal to 2 the the current
        %edge is a DOF and must be numbered. we must ,however, take care to
        %not renumber the edge so we compare the triangle_index to the two
        %common_triangles, if the triangle_index is the same as the lower
        %common triangle, then the edge in question has not yet been
        %numbered and a number is assigned to it. The current direction is
        %then also set to +1 which is out of the triangle, if it has been
        %numbered the we must go and find the numberin the lower triangle i
        %in question and attach the number here, and the current direction
        %must be set to -1. Coming into the triangle
        
        
        if size(common_quads,1) == 2 %This is a DOFF
            
            %common_triangles will always be sorted so we can compare our
            %tri_index to the last number in common_triangles
            if common_quads(2) ~= quad_index %hasn't been numbered 
                DOFF_NUM = DOFF_NUM + 1;
                N = N + 1;
                basis_supports = [basis_supports;common_quads'];
                %simply number the triangle and set the direction to +1
                
                %set directional flag to +1
                quadElements{quad_index}(direction_select(counter)) = 1; 
                %The edge is a DOFF so give it a number
                quadElements{quad_index}(doff_select(counter)) = DOFF_NUM;
                
                if (mod(quad_index,vertices) == 1) % First quad at contour point
                    quadElements{quad_index}(direction_select(2)) = -1;
                end
                if (mod(quad_index,vertices) == 0) % Last quad
                    quadElements{quad_index}(direction_select(4)) = 1; 
                end
%                 
            else
                %we set the  direction to -1 and fetch the DOFF number from
                %the lower triangle 
                
                %set directional flag to -1
                quadElements{quad_index}(direction_select(counter)) = -1; 
                %The edge is a DOFF so fetch the DOFF number from the lower
                %numbered triangle, this is not trivial as the DOFF
                %triangles{tri_index}(doff_select(counter)) = triangles{common_triangles(1)}(doff_select(counter));
                
                %retrieve the lower triangle our selected edge is a part
                %of.
                previous_quad = quadElements{common_quads(1)}(1:4);
              
                %Now we must find the edge in the previous quad and
                %retrieve the DOF number
                if curr_edge == previous_quad([1 2])
                    % DOFF 1 (12) is then edge AB
                    quadElements{quad_index}(doff_select(counter)) = quadElements{common_quads(1)}(doff_select(1));
                                    elseif curr_edge == previous_quad([1 4])
%                 elseif curr_edge == previous_quad([2 3])
                    % DOFF 2 (11) is then edge BC
                    quadElements{quad_index}(doff_select(counter)) = quadElements{common_quads(1)}(doff_select(2));
                elseif curr_edge == previous_quad([4 3])
                    % DOFF 3 (10) is then edge CD
                    quadElements{quad_index}(doff_select(counter)) = quadElements{common_quads(1)}(doff_select(3));
%                 elseif curr_edge == previous_quad([1 4])
                                    elseif curr_edge == previous_quad([2 3])
                    % DOFF 4 (9) is then edge AD
                    quadElements{quad_index}(doff_select(counter)) = quadElements{common_quads(1)}(doff_select(4));
                end
                %
                if (mod(quad_index,vertices) == 1) % First quad at contour point
                    quadElements{quad_index}(direction_select(2)) = -1;
                end
                if (mod(quad_index,vertices) == 0)
                    quadElements{quad_index}(direction_select(4)) = 1; % Last quad
                end

            end
            
        else
            % If we are in here it means the edge is not a DOFF and we
            % can simply put a -1 in the numbering system and a +1 in the
            % direction system.
            
            %set directional flag to +1
            quadElements{quad_index}(direction_select(counter)) = 1; 
            %The edge is not a DOFF so set number to -1
            quadElements{quad_index}(doff_select(counter)) = -1; 
        end
    end
end
% quadElements{1}(direction_select(4)) = -1;
%NEED TO CHECK HOW EXPENSIVE THIS IS COMPUTATIONALLY
quadElements = cell2mat(quadElements);
