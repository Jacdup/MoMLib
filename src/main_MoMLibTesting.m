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
mesh_create_option       = 5;              % 1 : Square plate
                                           % 2 : NASTRAN
                                           % 3 : Junction mesh
                                           % 4 : Cylinder
                                           % 5 : Sphere
TextOn                   = true;           % mesh visualisation text
flag_mesh_refine_uniform = false;
h_split_num              = 3;
quad                     = 1;
MBF                      = 0;


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
    %     [node_coords,tri_nodes] = ImportTriangleMeshNastran('bowtie_1.nas');
    [node_coords, tri_nodes] = ImportTriangleMeshNastran('C:\Users\19083688\Desktop\Masters\FEKO Models\Meshes\Hollow_cylinder.nas');
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

% Mesh a hollow cylinder
if mesh_create_option == 4
    
    rho = 0.5;
    vertices = 12;
    % vertices =10 + (9); % Number of vertices
    %  Contour = [0 0 0; 1 0 0; 2 0 0; 3 0 0; 4 0 0.1; 5 0 0.2; 6 0.1 0.3; 7 0.2 0.4; 8 0.25 0.5; 9 0.2 0.6; 10 0.2 0.5; 11 0.2 0.4 ]; % Hardcode some contour
    %  Contour = [0 0 0; 1 0 0; 2 0 0; 3 0 0; 4 0 3];
    %  Contour = [0 0 0; 0.1 0 0; 0.2 0 0; 0.3 0 0; 0.4 0 0; 0.5 0 0];
    Contour = [0 0 0; 1 0 0; 2 0 0];
    Contour = RefineMesh(Contour,1);
    [node_coords,quad_nodes, elements] = QuadMesh_v5(Contour,vertices,rho);
    [tri_nodes,triangles] = QuadtoTri(elements);
%     PlotMesh(node_coords,elements, triangles)
end

% Mesh a sphere
if mesh_create_option == 5
    [node_coords, quad_elems] = MeshCube(10);
    quad_nodes = {};
    for node = 1:length(quad_elems)
      quad_nodes = [quad_nodes;quad_elems(node,1:4),[0 0 0 0 0 0 0 0]]; 
    end
end

% Visualize the raw mesh data:
%PlotTriangleMeshRaw(node_coords,tri_nodes,TextOn);

%---------------------------------------------------------------
% Mesh pre-processing and visualisation --- the end result here is <mesh_data>
% which contains all the necessary derived information from the mesh,
% required for assigning dofs.
%---------------------------------------------------------------

% Create the mesh data, which includes edge definitions and connectivity:
if quad == 1
    [quad_blah,N,basis_supports_quad] = EdgeCalcQuad(quad_nodes,1);
    obs_basis_select = {1:N};
    src_basis_select      = {1:N};
    [common_basis_functions, num_obs, num_src] = common_basis(obs_basis_select,src_basis_select);
    quad_dof_idx = (1:N)';
    [new_quads,new_quad_points, new_quad_N, quad_observer_map, quad_source_map] = QuadBasisFunctionSelect(quad_blah ,common_basis_functions,basis_supports_quad,node_coords);
    MODE = 1;
    mex GCC='/usr/bin/gcc-7' -R2018a src\\Jacques\qmom.c
    [Z_mat] = qmom(new_quad_points, new_quads, new_quad_N, FREQUENCY,int32(quad_observer_map),int32(quad_source_map), num_obs, num_src, MODE);
    V_vec = qfillPlane(new_quad_points, new_quads, new_quad_N, FREQUENCY,int32(quad_observer_map), num_obs,0, MODE);
    if MBF == 1
        [U_Mat] = SelectDOFMBF(basis_supports_quad, vertices);
        Z_mat_redu = U_Mat.'*Z_mat*U_Mat;
        V_vec_redu = U_Mat.'*V_vec;
        Alpha_Vec_Redu = Z_mat_redu\V_vec_redu;
        I_vec = U_Mat*Alpha_Vec_Redu;
    end
else
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
%         V_vec = qfillPlane(new_quad_points, new_quads, new_quad_N, FREQUENCY,int32(quad_observer_map), num_obs,incident_dir, MODE);
            V_vec = CalcExciteVecPlaneWave(reduced_tri_dofs,reduced_node_coords, observer_map, k0, E_scalfac, theta_inc, phi_inc, eta_pol);
    end
    
    % Fill the impedance matrix:
    % mex  GCC='/usr/bin/gcc-7' -R2018a src\mom.c
    MODE = 1; % set the quadrature accuracy
    tic
    %   [Z_MoM] = qmom(new_quad_points, new_quads, new_quad_N, FREQUENCY,int32(quad_observer_map),int32(quad_source_map), num_obs, num_src, MODE);
    
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
end
% Solve:
if MBF == 0
    I_vec = Z_mat\V_vec;
end

%---------------------------------------------------------------
% Calculate and plot the surface current distribution.
%---------------------------------------------------------------

% Calculate the surface current density at the three vertices of all mesh triangle:
if quad == 1
    num_quad = size(new_quads,1);
    [quads_vertices_currents] = CalcElementsCurrentsQuad(I_vec, new_quads, quad_dof_idx, num_quad, node_coords);
    
     % Calculate centroid and vertex current magnitudes:
    quad_currents_cent = zeros(num_quad,1);
    quad_currents_dir = zeros(num_quad,3);
    quad_currents_vert = zeros(num_quad,4);
    for ii = 1:num_quad
        J_centroid              = (1/4)*[1 1 1 1]*squeeze(quads_vertices_currents(ii,1:4,1:3)); % Add each component at all 3 vertices and average
        for kk = 1:4
             quad_currents_dir(ii,1) =  quad_currents_dir(ii,1) + quads_vertices_currents(ii,kk,1);
             quad_currents_dir(ii,2) =  quad_currents_dir(ii,2) + quads_vertices_currents(ii,kk,2);
             quad_currents_dir(ii,3) =  quad_currents_dir(ii,3) + quads_vertices_currents(ii,kk,3);
        end

        quad_currents_cent(ii,1) = sqrt(abs(J_centroid*J_centroid'));
        for jj = 1:4
            J_vertex                 = squeeze(quads_vertices_currents(ii,jj,1:3)).'; % non-conjugate transpose to make this a row vector
            quad_currents_vert(ii,jj) = sqrt(abs(J_vertex*J_vertex')); % Magnitude
        end
    end
    quad_currents_dir_real = real(quad_currents_dir);
    quad_currents_dir_imag = imag(quad_currents_dir);
    
else
    num_tri = size(mesh_data.tri_nodes,1);
    [triangles_vertices_currents] = CalcElementsCurrents(I_vec, dof_data.tri_dofs, dof_data.tri_dofs_idx, num_tri, mesh_data.node_coords);
    
    % Calculate centroid and vertex current magnitudes:
    tri_currents_cent = zeros(num_tri,1);
    tri_currents_dir = zeros(num_tri,3);
    tri_currents_vert = zeros(num_tri,3);
    for ii = 1:num_tri
        J_centroid              = (1/3)*[1 1 1]*squeeze(triangles_vertices_currents(ii,1:3,1:3)); % Add each component at all 3 vertices and average
        tri_currents_cent(ii,1) = sqrt(abs(J_centroid*J_centroid'));
         for kk = 1:3
             tri_currents_dir(ii,1) =  tri_currents_dir(ii,1) + triangles_vertices_currents(ii,kk,1);
             tri_currents_dir(ii,2) =  tri_currents_dir(ii,2) + triangles_vertices_currents(ii,kk,2);
             tri_currents_dir(ii,3) =  tri_currents_dir(ii,3) + triangles_vertices_currents(ii,kk,3);
         end
         tri_currents_dir = abs(tri_currents_dir);
        
        
        for jj = 1:3
            J_vertex                 = squeeze(triangles_vertices_currents(ii,jj,1:3)).'; % non-conjugate transpose to make this a row vector
            tri_currents_vert(ii,jj) = sqrt(abs(J_vertex*J_vertex')); % Magnitude
        end
    end
end


% Visualize the current:
% PlotCurrent3D(0,false,mesh_data.tri_nodes,mesh_data.node_coords,tri_currents_cent);
if quad == 1
    PlotCurrent3DQuad(1,false,quad_blah,node_coords,quad_currents_vert);
    PlotCurrentDir3D(0,quad_blah,node_coords,quad_currents_dir_real,quad_currents_dir_imag);
else
    PlotCurrent3D(1,false,mesh_data.tri_nodes,mesh_data.node_coords,tri_currents_vert);
    PlotCurrentDir3D(1,mesh_data.tri_nodes,mesh_data.node_coords,tri_currents_dir);
end

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
plane               = 'XZ'; % set the observation plane: 'XY' or 'XZ' or 'YZ'
dalta_angle_degrees = 2;    % degree increment for observation directions
if quad == 1
    [farfield_XY] = qfarfield(new_quad_points, new_quads, new_quad_N, FREQUENCY, I_vec, dalta_angle_degrees, plane);
else
    [farfield_XY] = farfield(reduced_node_coords, reduced_tri_dofs, reduced_num_dofs, FREQUENCY, I_vec, dalta_angle_degrees, plane);
end

disp('<farfield> routine still needs some generalisation work');
figure;
plot(farfield_XY(:,1),sqrt(farfield_XY(:,3).^2 + farfield_XY(:,5).^2));
% [E_field_FEKO] = feko_farfield_extract('C:\Users\19083688\Desktop\MoM Codes\V_3.6_feat_speed\Jacques\FEKO\Hollow_Cylinder_0.1.txt');
% [E_field_FEKO] = feko_farfield_extract('C:\Users\19083688\Desktop\MoM Codes\V_3.6_feat_speed\Jacques\FEKO\Hollow_Cylinder_ZPlaneWave.txt');
% [E_field_FEKO] = feko_farfield_extract('C:\Users\19083688\Desktop\Masters\MoMLib\src\Jacques\FEKO\hollow_cyl_endcaps_xz.txt');
% [E_field_FEKO] = feko_farfield_extract('C:\Users\19083688\Desktop\Masters\MoMLib\src\Jacques\FEKO\hollow_cyl_8m_xz.txt');
[E_field_FEKO] = feko_farfield_extract('C:\Users\19083688\Desktop\Masters\MoMLib\src\Jacques\FEKO\sphere_0.5.txt');
hold on
plot(E_field_FEKO(:,1)*(pi/180),sqrt(E_field_FEKO(:,3).^2 + E_field_FEKO(:,5).^2));
% 
%---------------------------------------------------------------
return;
%---------------------------------------------------------------
