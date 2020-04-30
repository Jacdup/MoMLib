function [] = PlotTriangleMeshDofsMBF(mesh_data,dof_data,TextOn, U_Mat)
% Plot a mesh of triangles with node and triangle and global dof numbers.
% mesh_data    : output of the function <CreateMeshData.m>
% dof_data     : output of the function <CreateBasisFunctions.m>
%
% 2019-12-16: Created. MMB.

% Plot the nodes and triangles:
% PlotTriangleMeshRaw(mesh_data.node_coords,mesh_data.tri_nodes,TextOn);

% Get edge midpoints:
numedges       = size(mesh_data.edges,1);
num_dofs       = size(dof_data.dofs_to_edges,1);
ALPHA          = 0.5;
edge_midpoints = ALPHA*mesh_data.node_coords(mesh_data.edges(:,1),:) + (1-ALPHA)*mesh_data.node_coords(mesh_data.edges(:,2),:);
jj = 0;
if TextOn
%     for col = 1:3:length(U_Mat(1,:))
        for ii = 1:num_dofs
            
            if U_Mat(ii,4) ~= 0
                thisedge = dof_data.dofs_to_edges(ii,1);
                jj = jj + 1;
                %              plot3(edge_midpoints(thisedge,1),edge_midpoints(thisedge,2),edge_midpoints(thisedge,3),'rs');
                text(edge_midpoints(thisedge,1),edge_midpoints(thisedge,2),edge_midpoints(thisedge,3),[' ',num2str(U_Mat(ii,4))],'HorizontalAlignment','left','Color','r');
                
            end
            
            
        end
%     end
end
