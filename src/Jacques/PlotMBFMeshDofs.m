function [] = PlotMBFMeshDofs(mesh_data,dof_data,TextOn, DOF_mat)
% Plot a mesh of triangles with node and triangle and global dof numbers.
% mesh_data    : output of the function <CreateMeshData.m>
% dof_data     : output of the function <CreateBasisFunctions.m>
%
% 2019-12-16: Created. MMB. 
eps = 0.1;

% Plot the nodes and triangles:
% PlotTriangleMeshRaw(mesh_data.node_coords,mesh_data.tri_nodes,false);
% figure
hold on
% Get edge midpoints:
numedges       = size(mesh_data.edges,1);
num_dofs       = size(dof_data.dofs_to_edges,1);
ALPHA          = 0.5;
edge_midpoints = ALPHA*mesh_data.node_coords(mesh_data.edges(:,1),:) + (1-ALPHA)*mesh_data.node_coords(mesh_data.edges(:,2),:);

if TextOn 
    for ii = 1:length(nonzeros(DOF_mat))
        thisedge = dof_data.dofs_to_edges(DOF_mat(ii),1);
%         if (abs(edge_midpoints(thisedge,3) - 0.1) < eps) || (abs(edge_midpoints(thisedge,3)- 2.9375) < eps)
            plot3(edge_midpoints(thisedge,1),edge_midpoints(thisedge,2),edge_midpoints(thisedge,3),'rs');
            text(edge_midpoints(thisedge,1),edge_midpoints(thisedge,2),edge_midpoints(thisedge,3),[' ',num2str(DOF_mat(ii))],'HorizontalAlignment','left','Color','r');
%         end
    end
end
