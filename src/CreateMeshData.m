function [mesh_data] = CreateMeshData(node_coords,tri_nodes)
% This function creates a data structure containing all relevant
% information about the mesh, such that dofs can easily be assigned. 
%
% The results are returned in <mesh_data> (structure array data type) with 
% the following fields:
%
% mesh_data.node_coords          : node coordinates for each node 
%                                  ("points" in MATLAB terminology)
% mesh_data.tri_nodes            : triangle definitions with sorted nodes 
%                                  ("Triangulation connectivity list" in MATLAB terminology)
% mesh_data.connectivity_nodes   : Cell array of size num nodes; each cell
%                                  lists the triangles connected to that
%                                  node
% mesh_data.edges                : Array of size num global edges; each row
%                                  contais the two nodes making up a global
%                                  edge
% mesh_data.connectivity_edges   : Cell array of size num global edges; each
%                                  cell lists the triangles sharing that
%                                  edge
% mesh_data.tri_edges            : Array of size num tri; each row contains
%                                  the three, ordered global edge number,
%                                  with the local edges being defined as per
%                                  the local 3x2 matrix <local_edge_nodes_def>.  
%
% 2019-12-12: Created. MMB.

% Populate mesh_data.node_coords:
mesh_data.node_coords = node_coords;

% Populate mesh_data.tri_nodes:
mesh_data.tri_nodes = sort(tri_nodes,2);

% Create triangulation (MATLAB triangulation object of the mesh):
mesh_TR = triangulation(mesh_data.tri_nodes,mesh_data.node_coords);

% Populate mesh_data.connectivity_nodes:
mesh_data.connectivity_nodes = vertexAttachments(mesh_TR);
for ii = 1:size(mesh_data.connectivity_nodes,1)
    mesh_data.connectivity_nodes{ii} = sort(mesh_data.connectivity_nodes{ii},2);
end

% Populate mesh_data.edges:
mesh_data.edges = sort(edges(mesh_TR),2);

% Populate mesh_data.connectivity_edges:
mesh_data.connectivity_edges = edgeAttachments(mesh_TR,mesh_data.edges);
for ii = 1:size(mesh_data.connectivity_edges,1)
    mesh_data.connectivity_edges{ii} = sort(mesh_data.connectivity_edges{ii},2);
end

% Populate mesh_data.tri_edges:
numedges            = size(mesh_data.connectivity_edges,1);
numtri              = size(mesh_data.tri_nodes,1);
mesh_data.tri_edges = zeros(numtri,3);
for ii = 1:numedges
    conn_tri   = mesh_data.connectivity_edges{ii}; % all elements which share this edge
    numconn     = size(conn_tri,2);
    gledgenodes = mesh_data.edges(ii,:); % pair of ordered global node numbers of this edge
    for jj = 1:numconn % cycle over all connected elements and insert this edge number into <mesh_data.tri_edges>
        thistri                                = conn_tri(1,jj);
        thistrinodes                           = mesh_data.tri_nodes(thistri,:);
        localedge                              = LocalEdgeNum(gledgenodes,thistrinodes);
        mesh_data.tri_edges(thistri,localedge) = ii;
    end
end

%----------------------Internal functions----------------------
function [local_edge_num] = LocalEdgeNum(edge_nodes,tri_nodes)
% The local edge number is returned. Is is assumed that the input arguments
% are sorted lists of global node numbers, of lengths 2 and 3, respectively.
%
% 2019-12-13: Created. MMB.

%edge_select        = [1 2 ;1 3 ;2 3];
%direction_select   = [6 5 4];
%doff_select        = [9 8 7];

local_edge_nodes_def = [2 3
                        1 3
                        1 2];
if edge_nodes == tri_nodes(local_edge_nodes_def(1,:)) 
    local_edge_num = 1;
    return;
elseif edge_nodes == tri_nodes(local_edge_nodes_def(2,:)) 
    local_edge_num = 2;
    return;
else    
    local_edge_num = 3;
end

