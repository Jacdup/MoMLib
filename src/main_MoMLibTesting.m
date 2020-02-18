function [] = main_MoMLibTesting()
% This is the main file for testing and further development of the MoMLib
% package/toolbox. It is based upon the starting version established by
% Robey Beswick as part of his Master's thesis, Stellenbosch University, 2019.
%
% 2019-11-26: Created. MMB.

%---------------------------------------------------------------
% Initialize.
%---------------------------------------------------------------

%close all;

FREQUENCY                = 300e6;          % Solution frequancy
c0                       = 2.99792458e8;   % Speed of light in free space
lambda0                  = c0/FREQUENCY;   % Wavelength in free space
k0                       = 2*pi/lambda0;   % free space wavenumber
E_scalfac                = -1;              % plane wave spec
theta_inc                = 0;              % plane wave spec
phi_inc                  = 0;              % plane wave spec
eta_pol                  = pi/4;              % plane wave spec
flag_planewave           = true;
flag_lumped              = false;           % lumped sources and/or loads are present
interelem_VsrcZload      = [];             % init lumped circuit element specification
mesh_create_option       = 2;
TextOn                   = true;           % mesh visualisation text
flag_mesh_refine_uniform = false;
h_split_num              = 3;


%oldpath = path;
%path(oldpath,'..\V_6_feat_MemoryFix');

%---------------------------------------------------------------
% Read/create the mesh and visualise --- the end result here is <node_coords> 
% and <tri_nodes>, wherever they were obtained from.
%---------------------------------------------------------------

% Create a plate using basic MATLAB commands:
if mesh_create_option == 1
    coords_x            = [0:0.05:3];
    coords_y            = [0:0.02:0.5];
    [coords_x,coords_y] = meshgrid(coords_x,coords_y);
    node_coords         = [coords_x(:) coords_y(:)]; % convert matrices to column vectors
    node_coords(:,3)    = 0; % z-coord
    tri_nodes           = delaunay(node_coords(:,1:2));
end

% Import from a Nastran file that was created by FEKO:
if mesh_create_option == 2
    [node_coords,tri_nodes] = ImportTriangleMeshNastran('bowtie_1.nas');
    %node_coords(:,3)        = 0.1*rand(size(node_coords,1),1); % just a test
    interelem_VsrcZload     = [ 43 108 1 0 50 0]; % [tri- tri+ V_src^real V_src^imag Z_load^real Z_load^imag]
    flag_lumped             = true;
%     interelem_VsrcZload     = [ 43 108 0 0  5e4 0
%                                 90  80 1 0   0 0 
%                                 82  79 1 0   0 0]; % [tri- tri+ V_src^real V_src^imag Z_load^real Z_load^imag]
end

% Create a custom, junction mesh:
if mesh_create_option == 3
    node_coords = [0   0   0
                   1   0   0
                   1   1   0
                   0   1   0
                   0.5 0.5 1];
    tri_nodes   = [1   2   3
                   1   3   4
                   5   1   3];
end

% Visualize the raw mesh data:
%PlotTriangleMeshRaw(node_coords,tri_nodes,TextOn);

%---------------------------------------------------------------
% Mesh pre-processing and visualisation --- the end result here is <mesh_data> 
% which contains all the necessary derived information from the mesh,
% required for assigning dofs.
%---------------------------------------------------------------

% Create the mesh data, which includes edge definitions and connectivity:
[mesh_data] = CreateMeshData(node_coords,tri_nodes);
clear node_coords tri_nodes;

% Visualize the processed mesh data:
%PlotTriangleMeshProcessed(mesh_data,TextOn);

%---------------------------------------------------------------
% Make uniform mesh refinement.
%---------------------------------------------------------------

if flag_mesh_refine_uniform
    [node_coords,tri_nodes] = MeshRefineUniformSplit(mesh_data,h_split_num);
    [mesh_data] = CreateMeshData(node_coords,tri_nodes);
    clear node_coords tri_nodes;
    PlotTriangleMeshProcessed(mesh_data,TextOn);
end

%---------------------------------------------------------------
% Define all the RWG basis functions of the mesh.
%---------------------------------------------------------------

% Assign the dofs and define each associated basis function:
[dof_data,num_dofs] = CreateBasisFunctions(mesh_data);
%PlotTriangleMeshDofs(mesh_data,dof_data,TextOn);

%---------------------------------------------------------------
% Reduce the matrix if required; set up the excitation vector, system
% matrix, and solve. 
%---------------------------------------------------------------

% Build the sparse contributions to the full system due to lumped elements.
% The contributions can then be suitably reduced according to the selected
% sources and observers. The reason for not filling this after reduction,
% is because the lumped elements are specified in terms of the original
% mesh. Therefore the other option would be to recreate the lumped element
% specs in terms of the reduced mesh first and then work with that, but
% such an approach is more complicated, since it would effectively require
% that the full <dof_data> be reduced, since that is what is needed to
% assemple the lumped contributions (they are defined in terms of triangle
% numbers, not just i.t.o. dof numbers and triangle coords).
if flag_lumped
    [Zlump_rowcolval, Vlump_rowcolval] = CalcZmatVvecLumped(dof_data, num_dofs, mesh_data, interelem_VsrcZload);
end

% Code to reduce the number of dofs and also to selectively calculate matrix entries/blocks.
% Define observers and sources required:
% (V 3.1 allows for basis function selection, the result yields a N_new by
% N_new matrix, every basis function is both source and observer.)
obs_basis_select = [1:num_dofs];
src_basis_select = [1:num_dofs];

% Find the union of the observer and source list and append flags + record
% total observers and sources to calculate: 
[common_basis_functions, num_obs, num_src, global_to_redu_dofs] = CommonBasis(obs_basis_select,src_basis_select,num_dofs,true);

% Set up reduced set up dofs and mesh nodes, only relevant to the requested
% common basis functions:
[reduced_tri_dofs,reduced_node_coords, reduced_num_dofs, ...
 observer_map, source_map] = BasisFunctionSelect(dof_data.tri_dofs,common_basis_functions,dof_data.basis_supports,mesh_data.node_coords);

% Fill the excitation vector:
V_vec     = zeros(num_dofs,1);
if flag_planewave 
    V_vec = CalcExciteVecPlaneWave(reduced_tri_dofs,reduced_node_coords, observer_map, k0, E_scalfac, theta_inc, phi_inc, eta_pol);
end

% Fill the impedance matrix:
MODE = 1; % set the quadrature accuracy
tic
[Z_mat] = mom(reduced_node_coords, reduced_tri_dofs, reduced_num_dofs, FREQUENCY,int32(observer_map),int32(source_map), num_obs, num_src, MODE);
toc

% Add lumped circuit element contributions to the reduced Z and V, by
% obtaining the contributions in sparse format and then adding to the
% system:  
if flag_lumped
    [Zlump_reduced, Vlump_reduced] = FinalZmatVvecLumped(Zlump_rowcolval, Vlump_rowcolval, global_to_redu_dofs, num_obs, num_src, observer_map, source_map);
    Z_mat = Z_mat + Zlump_reduced;
    V_vec = V_vec + Vlump_reduced;
end

% Solve:
I_vec = Z_mat\V_vec;

%---------------------------------------------------------------
% Calculate and plot the surface current distribution.
%---------------------------------------------------------------

% Calculate the surface current density at the three vertices of all mesh triangle:
num_tri = size(mesh_data.tri_nodes,1);
[triangles_vertices_currents] = CalcElementsCurrents(I_vec, dof_data.tri_dofs, dof_data.tri_dofs_idx, num_tri, mesh_data.node_coords);

% Calculate centroid and vertex current magnitudes:
tri_currents_cent = zeros(num_tri,1);
tri_currents_vert = zeros(num_tri,3);
for ii = 1:num_tri
    J_centroid              = (1/3)*[1 1 1]*squeeze(triangles_vertices_currents(ii,1:3,1:3));
    tri_currents_cent(ii,1) = sqrt(abs(J_centroid*J_centroid'));
    for jj = 1:3
        J_vertex                 = squeeze(triangles_vertices_currents(ii,jj,1:3)).'; % non-conjugate transpose to make this a row vector
        tri_currents_vert(ii,jj) = sqrt(abs(J_vertex*J_vertex'));
    end
end

% Visualize the current:
%PlotCurrent3D(0,false,mesh_data.tri_nodes,mesh_data.node_coords,tri_currents_cent);
PlotCurrent3D(1,false,mesh_data.tri_nodes,mesh_data.node_coords,tri_currents_vert);

%---------------------------------------------------------------
% Custom impedance calculation for the 'bowtie_1.nas' mesh, when excited at the neck: edge nodes 24 and 28; dof number 49.
%---------------------------------------------------------------

if flag_lumped
    if num_dofs > 48
        edgevectemp = mesh_data.node_coords(28,:) - mesh_data.node_coords(24,:);
        ELength     = sqrt(edgevectemp*edgevectemp');
        Z_11        = 1/(ELength*I_vec(49,1))
        disp('General source impedance calculation must still be implemented');
    end
end

%---------------------------------------------------------------
% Calculate and plot the farfield.
%---------------------------------------------------------------

% The output of the farfield routine has this row format:
% THETA     PHI     magn{Etheta} phase{Etheta} magn{Ephi} phase{Ephi}
%   1        2            3           4            5          6
plane               = 'XY'; % set the observation plane: 'XY' or 'XZ' or 'YZ'
dalta_angle_degrees = 2;    % degree increment for observation directions
[farfield_XY] = farfield(reduced_node_coords, reduced_tri_dofs, reduced_num_dofs, FREQUENCY, I_vec, dalta_angle_degrees, plane);
disp('<farfield> routine still needs some generalisation work');
figure;
semilogy(farfield_XY(:,2),sqrt(farfield_XY(:,3).^2 + farfield_XY(:,5).^2));

%---------------------------------------------------------------
return;
%---------------------------------------------------------------
