function [triangles,N, basis_supports] = EdgeCalc(triangles, order)%filename recieved is a string
%Robey Beswick
%August 2018
%18472648@sun.ac.za
%
%DESCRIPTION
%this routine allocates the current directions over the entire mesh. It
%must be done flawlessly every time. The current is returned using the
%current direction is defined over the edge as A B C a b c, where 'a'
%determines the current direction over the BC '1' means the current flows
%out and '-1' means the current flows in.
%
%This functions allocates the current directions and DOFF numbering over
%the triangular mesh. It recieves a sorted triangles cell matrix.
%The layout is :
%       A  B  C || D  E  F || H  I  J
%       
%       current directions:
%       (+1) current going out (-1) coming in
%       D -> current direction over edge BC
%       E -> current direction over edge AC
%       F -> current direction over edge AB
%
%       DOF numbering:
%       if (-1) not a DOF otherwise DOF number
%       H -> number for edge BC
%       I -> number for edge AC
%       J -> number for edge AB
%
%UPDATE 24/04/2019
%The code is updated from version 3.0 to now also return the basis function
%supports, an N by 2 array that contains the two triangles supporting the
%basis function at the nth DOFF(the row index is the doff number it supports)
%This is a trivial implementation.

%Define variables and constants
connectivity_data = {};
edge_select = [1 2 ;1 3 ;2 3];
direction_select = [6 5 4];
doff_select = [9 8 7];
local_edge_select = [1,2];
DOFF_NUM = 0;
tri = [0, 0, 0];
high_num = 0;
N = 0;
basis_supports = [];
add_tri = 0;

%populate connectivity list data
for tri_index = 1:size(triangles,1)
    tri = triangles{tri_index,1:3};
    for node = 1:3
      
        if high_num >= tri(node)%means we dont have to add another cell to the end
            
            connectivity_data{tri(node)} = [connectivity_data{tri(node)};tri_index];
        else %new node, create new cell
            
            connectivity_data{tri(node)} = [tri_index];
            high_num = tri(node);
        end
    end
end
% connectivity_data{1}
% connectivity_data{2}
% connectivity_data{3}
% connectivity_data{4}
% connectivity_data{5}
% connectivity_data{6}


%-----------%-----------%-----------%-----------%-----------%-----------%-----------%-----------
%-----------Now we will simultaneously number the DOF's and assign the current-----------
%-----------directions for the triangles.-----------
%-----------%-----------%-----------%-----------%-----------%-----------%-----------%-----------
% loc_num = [];
% for tri_index = 1:(size(triangles,1))
%      tri = triangles{tri_index,1:3};
% %     for node = 1:3
%    
%         curr_edge(tri_index,1:2) = [tri(edge_select(1,1)),tri(edge_select(1,2))];
%         curr_edge(tri_index,3:4) = [tri(edge_select(2,1)),tri(edge_select(2,2))];
%         curr_edge(tri_index,5:6) = [tri(edge_select(3,1)),tri(edge_select(3,2))];
% %          loc_num(tri_index,1:6) = [curr_edge; [1 ,2,1,2,1,2]];
%         [tri(1) tri(2) tri(3)]
%         [1 2; 1 2; 1 2]
%         
%         tri_struct1:
%         [tri(1) tri(2) 1 2] (edge_select(1))
%         [tri(2) tri(3) 1 2] (edge_select(3))
%         [tri{3) tri(1) 1 2] (edge_select(2)) -> NB: this is switched!
%         
%         tri_struct2:
%         [tri2(1) tri2(2) 1 2]
%         [...]
%         [...]
%         
%         curr_edge = [2 3] 
%         
%         if tri_struct1(edge_select(counter,1)) == tri_struct2(edge_select(counter,1)
%         % I just have to assign two numbers (1/2) to each local node
%         
%         
% %     end
% end


%Now we will iterate through the triangles 

for tri_index = 1:(size(triangles,1))
    tri = triangles{tri_index,1:3};
    for counter = 1:3 % Edges
        if order == 1 
            first_order = mod(tri_index, 2); % Every second element is first order
        else
            first_order = 1;
        end
        %once we have selected a triangle T_n then we will select an edge
        %from that triangles
        curr_edge = [tri(edge_select(counter,1)),tri(edge_select(counter,2))]; % Global edge
%         curr_edge_local = [curr_edge;local_edge_select]; % Local edge
        
        %once the edge of the triangle is selected then we go the the
        %connectivity_list and check if the edge shares two triangles if
        %the edges does indeed share two triangles then it is a DOF. see
        %below.
        
        %choose the relevant connectivity list columns
        node_A = connectivity_data{curr_edge(1)};
        node_B = connectivity_data{curr_edge(2)};
        common_triangles = intersect(node_A,node_B);
        [local_edge_num1] = LocalEdgeNum(curr_edge,triangles{common_triangles(1)}(1:3)); % This is the first triangles number
        if size(common_triangles,1) ==2
            ind = 2;
        else
            ind = 3;
        end
        [local_edge_num2] = LocalEdgeNum(curr_edge,triangles{common_triangles(ind)}(1:3)); % This is the second triangles number
        
            
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
        
        % Edit 29/01/2020
        if order == 0
            comm_size = 2;
        else
            comm_size = 4;
        end
        if size(common_triangles,1) == comm_size %This is a DOFF
            
            %common_triangles will always be sorted so we can compare our
            %tri_index to the last number in common_triangles
            if order == 1
               true_flag = (common_triangles(3) ~= tri_index && common_triangles(4) ~= tri_index);
            else
                true_flag = (common_triangles(2) ~= tri_index);
            end
            if  true_flag%hasn't been numbered 
%                 if first_order == 1
                    DOFF_NUM = DOFF_NUM + 1;
                    N = N + 1;
%                 end
                
                basis_supports = [basis_supports;common_triangles'];
                %simply number the triangle and set the direction to +1
                
                %set directional flag to +1
                triangles{tri_index}(direction_select(counter)) = 1;
                if first_order == 0 %Confusing, but this means it is first order.
                    if local_edge_num1(1) == local_edge_num2(1)
                        triangles{tri_index}(direction_select(counter)) = 2;
                    else
                    
%                                         if curr_edge
                    triangles{tri_index}(direction_select(counter)) = 2;
                    end
                    %                     DOFF_NUM = DOFF_NUM - 1;
                    %                 else
                    %                     DOFF_NUM = DOFF_NUM + 1;
                end
                    
                %The edge is a DOFF so give it a number
                triangles{tri_index}(doff_select(counter)) = DOFF_NUM;
%                 if first_order == 0 
%                     triangles{tri_index}(doff_select(counter)) = triangles{tri_index-1}(doff_select(counter));
% %                     triangles{tri_index}(doff_select(counter)) = DOFF_NUM -1 ; 
%                 end
                
            else % Already has a DOF
                %we set the  direction to -1 and fetch the DOFF number from
                %the lower triangle 
                
                %set directional flag to -1
                triangles{tri_index}(direction_select(counter)) = -1; 
                if first_order == 0 
                      if local_edge_num1(1) == local_edge_num2(1)
                        triangles{tri_index}(direction_select(counter)) = -2;
                      else
                    %                     if curr_edge
                        triangles{tri_index}(direction_select(counter)) = 2;
                      end
                end
                %The edge is a DOFF so fetch the DOFF number from the lower
                %numbered triangle, this is not trivial as the DOFF
                %triangles{tri_index}(doff_select(counter)) = triangles{common_triangles(1)}(doff_select(counter));
                
                %retrieve the lower triangle our selected edge is a part
                %of.
%                 for ii=1:2 
                  if first_order == 0
                      add_tri = 1;
                  else
                      add_tri = 0;
                  end
                    previous_triangle = triangles{common_triangles(1+add_tri)}(1:3);

                    %Now we must find the edge in the previous triangle and 
                    %retrieve the DOF number
                    if curr_edge == previous_triangle([1 2])
                        triangles{tri_index}(doff_select(counter)) = triangles{common_triangles(1+add_tri)}(doff_select(1));
                    elseif curr_edge == previous_triangle([1 3])
                        triangles{tri_index}(doff_select(counter)) = triangles{common_triangles(1+add_tri)}(doff_select(2));
                    elseif curr_edge == previous_triangle([2 3])
                        triangles{tri_index}(doff_select(counter)) = triangles{common_triangles(1+add_tri)}(doff_select(3));
                    end
%                   end
                    
            end
            
        else
            % If we are in here it means the edge is not a DOFF and we
            % can simply put a -1 in the numbering system and a +1 in the
            % direction system.
            
            %set directional flag to +1
            triangles{tri_index}(direction_select(counter)) = 1; 
            if first_order == 0 
                 if local_edge_num1 == local_edge_num2
                        triangles{tri_index}(direction_select(counter)) = 2;
                    else
                    
                    %                     if curr_edge
                    triangles{tri_index}(direction_select(counter)) = 2;
                  end
%                     triangles{tri_index}(direction_select(counter)) = 2; 
            end
            %The edge is not a DOFF so set number to -1
            
            triangles{tri_index}(doff_select(counter)) = -1; 
%             if first_order == 0 
%                     triangles{tri_index}(direction_select(counter)) = 2; 
%              end
        end
    end
end

%NEED TO CHECK HOW EXPENSIVE THIS IS COMPUTATIONALLY
triangles = cell2mat(triangles);
