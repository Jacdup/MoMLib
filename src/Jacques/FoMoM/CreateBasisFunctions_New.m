function [dof_data,num_dofs] = CreateBasisFunctions_New(mesh_data, first_order)
% This function assigns RWG degrees of freedom (dofs) to all shared edges of
% the entire mesh <tri_nodes> and assigns the current directions over each
% edge, such that the RWGs are uniquely defined. Unshared edges do not get
% dofs. The function supports meshes with junctions (edges shared by more
% than two triangles). In the special case of no junctions (e.g. a flat plate)
% the dimension of <tri_dofs> will be the same or less than num triangles. It
% can be less when a given triangles has no dofs. With junctions, its
% dimension can be more than num triangles, but always less than 2*numtri.
% At a junction a given triangle can have a maximum of two basis functions
% associated with its junction edge.
%
% The results are returned in <tri_dofs> and <basis_supports>.
%
% Definition of values in a row of <tri_dofs>:
%       A  B  C || D  E  F || G  H  I
%
%       triangle node numbers:
%       A -> node 1
%       B -> node 2
%       C -> node 3
%
%       current directions:
%       +1 ==> current flows out of this triangle; -1 ==> current flows into this triangle
%       (default initialization values are +1)
%       D -> current direction over edge BC
%       E -> current direction over edge AC
%       F -> current direction over edge AB
%
%       DOF numbering:
%       if (-1) not a DOF otherwise DOF number
%       G -> number for edge BC
%       H -> number for edge AC
%       I -> number for edge AB
%
% Definition of values in a row of <basis_supports>, one row for each RWG dof:
%       [  triangle_RWG_current_flows_out_of  triangle_RWG_current_flows_into  ]
%
% Definition of <tri_dofs_idx>. Column vector with length equal to size of
% <tri_dofs>. Each row contain one value:
% [  global_triangle_number_of_this_row_in_<tri_dofs>  ]
%
% Definition of <dofs_to_edges>. Column vector with length equal to size of
% <num_dofs>. Each row contain one value:
% [  global_edge_associated_with_this_dof  ]
%
% Definition of <tri_to_dofs>. Struct of size num tri, with the dofs
% associated with each tri lested in a struct element <dofs>.
%
% 2018-08-00: Created by Robey Beswick, 18472648@sun.ac.za.
% 24/04/2019: UPDATE. The code is updated from version 3.0 to now also
% return the basis function supports, an N by 2 array that contains the two
% triangles supporting the basis function at the nth DOFF(the row index is
% the doff number it supports) This is a trivial implementation. RCB.
% 2019-12-12: Renamed from <EdgeCalc>. Polished up and added node sorting and moved expansion of
% "triangles" with 6 zeros after the node data, to here, instead of having
% it as part of mesh importing routine. MMB.
% 2019-12-12: Rewrote completely, based on comprehensive connectivity data,
% such that dofs are assigned per edge. Each global edge is processed in
% turn. MMB.
% 2019-12-16: Added <tri_dofs_idx> and <dofs_to_edges> and encapsulated all
% the dof date in <dof_data>. MMB.

% Define variables and constants:
numtri             = size(mesh_data.tri_nodes,1);
numedg             = size(mesh_data.edges,1);

% Init the output:
tri_dofs       = [];
basis_supports = zeros(3*numtri,2); % preallocate for speed, unreachable upper estimate --> two dofs per shared edge, i.e. (numtri x 3 x 2) / 2
dofs_to_edges  = zeros(3*numtri,2);
num_dofs       = 0;

% Assign dofs and creat basis functions by filling a structure with the
% values and then converting the structure to a standard array containing
% only the relevant nonzero entries. The structure allows for assigning two
% dofs to each triangle edge, which is the absolute maximum possible.

% Update 2020/04/06: Maintaining max of 2 DOFs per edge, but now with first
% order. If 3 DOFs are required as at a junction, this still needs to be
% implemented.

% Structure of a struct entry (there are numtri entries):
% item 1: nnz (how many nonempty rows in the second item: 0,1 or 2)
% item 2: [D  E  F  G  H  I; D  E  F  G  H  I] which allows for associating two dofs with each edge of a triangle
for ii = numtri:-1:1 % change the order, so that the array is preallocated once
    tri_dofs_struct(ii).nnz  = 0;                                 % init the struct
    tri_dofs_struct(ii).dofs = [1 1 1 -1 -1 -1; 1 1 1 -1 -1 -1];  % init the struct
end

% Now populate <tri_dofs_struct> by cycling through all global edges:
for ii = 1:numedg
    edge_tris = mesh_data.connectivity_edges{ii}; % use curly brackets to get the actual contents of the cell (<connectivity_edges> is a cell array)
    num_conn  = size(edge_tris,2) - first_order; % number of triangles sharing this edge
    
    % Edit 2020/04/06
    %     edge_tris = [mesh_data.connectivity_edges{ii}, mesh_data.connectivity_edges{ii}(1)];
    %     num_conn  = size(edge_tris,2); % All triangles now have an extra DOF
    
    %     if num_conn > 1 % then one or more dofs will have to be assigned to this edge
    if num_conn > (1+first_order) % then one or more dofs will have to be assigned to this edge
        edge_nodes = mesh_data.edges(ii,:);
        for jj = 1:num_conn-1 % one less dof assigned to a junction, than the number of currents meeting there (KCL)
            
            num_dofs                     = num_dofs + 1;
            dofs_to_edges(num_dofs,1)    = ii;
            t1                           = edge_tris(jj);
            t2                           = edge_tris(jj+1+first_order);
            basis_supports(num_dofs,1:2) = [t1 t2];
            
            % Add dof data for t1 (current flows out):
            e1 = LocalEdgeNum(edge_nodes,mesh_data.tri_nodes(t1,:));
            d1 = LocalEdgeOrientation(edge_nodes,mesh_data.tri_nodes(t1,:));
            if     tri_dofs_struct(t1).dofs(1,3+e1) == -1 % then add the dof to the first line of the 2x6 datablock
                tri_dofs_struct(t1).nnz          = max(tri_dofs_struct(t1).nnz,1);
                tri_dofs_struct(t1).dofs(1,3+e1) = num_dofs;
                
                if (mod(t1,2) == 0) && first_order % Every second element in first order construction
                    tri_dofs_struct(t1).dofs(1,e1)   = 2; % not really necessary, as default is +1
                else
                    tri_dofs_struct(t1).dofs(1,e1)   = 1;
                end
                
            elseif tri_dofs_struct(t1).dofs(2,3+e1) == -1
                tri_dofs_struct(t1).nnz          = max(tri_dofs_struct(t1).nnz,2);
                tri_dofs_struct(t1).dofs(2,3+e1) = num_dofs;
                if (mod(t1,2) == 0) && first_order % Every second element in first order construction
                    tri_dofs_struct(t1).dofs(2,e1)   = 2; % not really necessary, as default is +1
                else
                    tri_dofs_struct(t1).dofs(2,e1)   = 1;
                end
                
            else
                %                 tri_dofs_struct(t1).dofs(2,3+e1) = -1; % Revert back, since this edge is probably shared
                error('Attempt at assigning a third dof to a triangle edge');
            end
            
            % Add dof data for t2 (current flows in):
            e2 = LocalEdgeNum(edge_nodes,mesh_data.tri_nodes(t2,:));
            d2 = LocalEdgeOrientation(edge_nodes,mesh_data.tri_nodes(t2,:));
            if     tri_dofs_struct(t2).dofs(1,3+e2) == -1 % then add the dof to the first line of the 2x6 datablock
                tri_dofs_struct(t2).nnz          = max(tri_dofs_struct(t2).nnz,1);
                tri_dofs_struct(t2).dofs(1,3+e2) = num_dofs;
                
                if (mod(t2,2) == 0) && first_order % Every second element in first order construction
                    if d1(1) ~= d2(1)
                        tri_dofs_struct(t2).dofs(1,e2)  = 2;
                    else
                        tri_dofs_struct(t2).dofs(1,e2)  = -2;
                    end
                else
                    tri_dofs_struct(t2).dofs(1,e2)   = -1;
                end
                
            elseif tri_dofs_struct(t2).dofs(2,3+e2) == -1
                tri_dofs_struct(t2).nnz          = max(tri_dofs_struct(t2).nnz,2);
                tri_dofs_struct(t2).dofs(2,3+e2) = num_dofs;
                if (mod(t2,2) == 0) && first_order
                    if d1(1) ~= d2(1)
                        tri_dofs_struct(t2).dofs(2,e2)   = 2;
                    else
                        tri_dofs_struct(t2).dofs(2,e2)  = -2;
                    end
                else
                    tri_dofs_struct(t2).dofs(2,e2)   = -1;
                end
            else
                error('Attempt at assigning a third dof to a triangle edge');
            end
            
        end
    end
end
first_order = 0;
% Assign final outputs:
basis_supports = basis_supports(1:num_dofs,:); % trim
dofs_to_edges  = dofs_to_edges(1:num_dofs,:);  % trim

% Compose the trimmed, standard 2D array data block:
count_doftris = 0;
for ii = 1:1+first_order:numtri
    count_doftris = count_doftris + tri_dofs_struct(ii).nnz;
end
tri_dofs     = zeros(count_doftris,9);
tri_dofs_idx = zeros(count_doftris,1);
count_doftris = 0;
for ii = 1+first_order:1+first_order:numtri
    for jj = 1:tri_dofs_struct(ii).nnz
        count_doftris                 = count_doftris + 1;
        tri_dofs(count_doftris,1:3)   = mesh_data.tri_nodes(ii,:);
        tri_dofs(count_doftris,4:9)   = tri_dofs_struct(ii).dofs(jj,:);
        tri_dofs_idx(count_doftris,1) = ii;
    end
end

% Create lookup data tri-to-dofs, to thsart the dofs associated with a
% given triangles is readily available. This is required when applying
% voltage sources and lumped loads, which are specified at pairs of
% triangles. I.e. this data is required to at O(1) cost, given a pair of
% triangles, find the corresponding dof (which is the row in
% basis_supports). A struct datastructure is used, as it allows for a
% element of varying size.
numpertri                = zeros(numtri); % init
tri_to_dofs(numtri).dofs = [];            % preallocate for speed
for ii = 1:num_dofs
    for jj = 1:2
        thistri                                       = basis_supports(ii,jj);
        numpertri(thistri)                            = numpertri(thistri) + 1;
        tri_to_dofs(thistri).dofs(numpertri(thistri)) = ii;
    end
end

% Assign encapsulated output data struct:
dof_data.tri_dofs       = tri_dofs;
dof_data.tri_dofs_idx   = tri_dofs_idx;
dof_data.dofs_to_edges  = dofs_to_edges;
dof_data.basis_supports = basis_supports;
dof_data.tri_to_dofs    = tri_to_dofs;

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

function [local_edge_dir] = LocalEdgeOrientation(edge_nodes,tri_nodes)
% The local edge number is returned. Is is assumed that the input arguments
% are sorted lists of global node numbers, of lengths 2 and 3, respectively.
%
% 2019-12-13: Created. MMB.

%edge_select        = [1 2 ;1 3 ;2 3];
%direction_select   = [6 5 4];
%doff_select        = [9 8 7];

local_edge_nodes_def = [1 2
    1 3
    2 3];
if edge_nodes == tri_nodes(local_edge_nodes_def(1,:))
    local_edge_dir = [1 2];
    return;
elseif edge_nodes == tri_nodes(local_edge_nodes_def(2,:))
    local_edge_dir = [2 1];
    return;
else
    local_edge_dir = [1 2];
end

