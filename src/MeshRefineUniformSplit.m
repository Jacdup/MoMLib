function [node_coords,tri_nodes] = MeshRefineUniformSplit(mesh_data,h_split_num)
% This function recieves a mesh with its connectivity data. From that a new
% mesh is created, which is the result of splitting each original triangle
% edge into <h_split_num> equal segments (with the necessary internal nodes
% added as dictated by the value of <h_split_num>). Each original triangle
% is then split into a set of identical, congruent, new mesh triangles,
% with a shared global nodes list such that shared nodes are correctly
% referenced on multiple instances by all the relevant triangles.    
%
% See MMB notes in the \doc folder for details to make sense of the
% implementation. Not necessary to understand in order to use this
% routine, of course.
%
% The input and output arguments are as follows:
%
% mesh_data     : this is the original mesh and its connectivity data, as established in routine <CreateMeshData.m> (structure array data type)
% h_split_num   : edge splitting factor --- 2 and 3 are supported
% node_coords   : node coordinates for each node in the new mesh (x,y,z coordinates in each row) 
% tri_nodes     : three nodes defining each triangle in the new mesh, is listed in each row (unsorted) 
%
% 2019-12-24: Created. MMB.

% Init:
num_tri_old   = size(mesh_data.tri_nodes,1);
num_edges_old = size(mesh_data.edges,1);
num_nodes_old = size(mesh_data.node_coords,1);
if h_split_num == 2 || h_split_num == 3 % if equal to 2 or to 3
    tri_nodes                      = zeros( (h_split_num^2)*num_tri_old , 3 );
    node_coords                    = zeros( num_nodes_old + (h_split_num-1)*num_edges_old + ((h_split_num-2)^2)*num_tri_old , 3 );
    node_coords(1:num_nodes_old,:) = mesh_data.node_coords; % keep the old nodes, new ones will be appended to the list
else
    error('Invalid <h_split_num> in MeshRefineUniformSplit');
end

% Create new edge and face nodes numbers and append their coordinates to
% <node_coords> as initialized above, to form the new <node_coords>: 
if h_split_num == 2

    % New edge nodes:
    new_edgenodes = num_nodes_old + [1:num_edges_old]'; % for each old global edge, each row contains 
                                                        % the ordered new global node number, subdividing 
                                                        % that edge uniformly   
    node_coords(num_nodes_old+1:num_nodes_old+num_edges_old ,:) = ...
        0.5*mesh_data.node_coords( mesh_data.edges(:,1) ,:) +  ...
        0.5*mesh_data.node_coords( mesh_data.edges(:,2) ,:);

    % New face nodes:
    % (none in the case of h_split_num = 2)

elseif h_split_num == 3
    
    % New edge nodes:
    new_edgenodes = num_nodes_old + [ [1:num_edges_old]' (num_edges_old + [1:num_edges_old]') ]; % two new nodes per edge
    node_coords(num_nodes_old+1:num_nodes_old+num_edges_old ,:) = ... % first edge nodes, closest to the first node edge
        (2/3)*mesh_data.node_coords( mesh_data.edges(:,1) ,:) +  ...
        (1/3)*mesh_data.node_coords( mesh_data.edges(:,2) ,:);
    node_coords(num_nodes_old+num_edges_old+1:num_nodes_old+2*num_edges_old ,:) = ... % second edge nodes, closest to the second node edge
        (1/3)*mesh_data.node_coords( mesh_data.edges(:,1) ,:) +  ...
        (2/3)*mesh_data.node_coords( mesh_data.edges(:,2) ,:);
    
    % New face nodes:
    new_facenodes = num_nodes_old + 2*num_edges_old + [1:num_tri_old]'; % one new node per face/tri/elem
    node_coords(num_nodes_old+2*num_edges_old+1:num_nodes_old+2*num_edges_old+num_tri_old ,:) = ... % centroids of the triangles
        (1/3)*mesh_data.node_coords( mesh_data.tri_nodes(:,1) ,:) +  ...
        (1/3)*mesh_data.node_coords( mesh_data.tri_nodes(:,2) ,:) +  ...
        (1/3)*mesh_data.node_coords( mesh_data.tri_nodes(:,3) ,:);

end

% Fill <tri_nodes> by splitting each old triangle into h_split_num^2 new
% ones. To make sense of this part, relies heavily on MMB notes 2019-12-24,
% available in the MoMLib\doc folder: 
if h_split_num == 2

    for ii = 1:num_tri_old
    
        % Create a map from a local node grid to the global nodes. The new
        % triangles are defined in terms of this local node grid:
        localgridnode( 1) = mesh_data.tri_nodes(ii,1);
        tempedge          = mesh_data.tri_edges(ii,3);
        localgridnode( 2) = new_edgenodes(tempedge,1);
        localgridnode( 3) = mesh_data.tri_nodes(ii,2);
        tempedge          = mesh_data.tri_edges(ii,2);
        localgridnode( 4) = new_edgenodes(tempedge,1);
        tempedge          = mesh_data.tri_edges(ii,1);
        localgridnode( 5) = new_edgenodes(tempedge,1);
        localgridnode( 6) = mesh_data.tri_nodes(ii,3);
        
        % Create the set of new triangles associated with this old element:
        tri_nodes( (ii-1)*(h_split_num^2)+1:ii*(h_split_num^2) , : ) = ... 
            [ localgridnode( 1)  localgridnode( 2)  localgridnode( 4)
              localgridnode( 2)  localgridnode( 3)  localgridnode( 5)
              localgridnode( 4)  localgridnode( 2)  localgridnode( 5)
              localgridnode( 4)  localgridnode( 5)  localgridnode( 6) ];
    end
    
elseif h_split_num == 3

    for ii = 1:num_tri_old
    
        % Create a map from a local node grid to the global nodes. The new
        % triangles are defined in terms of this local node grid:
        localgridnode( 1) = mesh_data.tri_nodes(ii,1);
        tempedge          = mesh_data.tri_edges(ii,3);
        localgridnode( 2) = new_edgenodes(tempedge,1);
        localgridnode( 3) = new_edgenodes(tempedge,2);
        localgridnode( 4) = mesh_data.tri_nodes(ii,2);
        tempedge          = mesh_data.tri_edges(ii,2);
        localgridnode( 5) = new_edgenodes(tempedge,1);
        localgridnode( 6) = new_facenodes(ii,1);
        tempedge          = mesh_data.tri_edges(ii,1);
        localgridnode( 7) = new_edgenodes(tempedge,1);
        tempedge          = mesh_data.tri_edges(ii,2);
        localgridnode( 8) = new_edgenodes(tempedge,2);
        tempedge          = mesh_data.tri_edges(ii,1);
        localgridnode( 9) = new_edgenodes(tempedge,2);
        localgridnode(10) = mesh_data.tri_nodes(ii,3);
        
        % Create the set of new triangles associated with this old element:
        tri_nodes( (ii-1)*(h_split_num^2)+1:ii*(h_split_num^2) , : ) = ... 
            [ localgridnode( 1)  localgridnode( 2)  localgridnode( 5)
              localgridnode( 2)  localgridnode( 3)  localgridnode( 6)
              localgridnode( 3)  localgridnode( 4)  localgridnode( 7)
              localgridnode( 5)  localgridnode( 2)  localgridnode( 6)
              localgridnode( 6)  localgridnode( 3)  localgridnode( 7)
              localgridnode( 5)  localgridnode( 6)  localgridnode( 8)
              localgridnode( 6)  localgridnode( 7)  localgridnode( 9)
              localgridnode( 8)  localgridnode( 6)  localgridnode( 9)
              localgridnode( 8)  localgridnode( 9)  localgridnode(10) ];
    end
end
