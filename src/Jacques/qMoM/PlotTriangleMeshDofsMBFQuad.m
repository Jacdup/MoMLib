function [] = PlotTriangleMeshDofsMBFQuad(node_coords,quad_blah,TextOn, U_Mat, edges, dofs_to_edges)
% Plot a mesh of triangles with node and triangle and global dof numbers.
% mesh_data    : output of the function <CreateMeshData.m>
% dof_data     : output of the function <CreateBasisFunctions.m>
%
% 2019-12-16: Created. MMB.

% Plot the nodes and triangles:
% PlotTriangleMeshRaw(mesh_data.node_coords,mesh_data.tri_nodes,TextOn);

% Get edge midpoints:
% numedges       = size(mesh_data.edges,1);
num_dofs       = size(dofs_to_edges,1);
% num_dofs       = max(max(quad_blah(:,9:12)));
ALPHA          = 0.5;

edge_midpoints = ALPHA*node_coords(edges(:,1),:) + (1-ALPHA)*node_coords(edges(:,2),:);
jj = 0;
if TextOn
    for col = [1,4]
        for ii = 1:num_dofs
            
            if U_Mat(ii,col) ~= 0
                thisedge = dofs_to_edges(ii,1);
                jj = jj + 1;
                %              plot3(edge_midpoints(thisedge,1),edge_midpoints(thisedge,2),edge_midpoints(thisedge,3),'rs');
                text(edge_midpoints(thisedge,1),edge_midpoints(thisedge,2),edge_midpoints(thisedge,3),[' ',num2str(U_Mat(ii,col+1))],'HorizontalAlignment','left','Color','r');
                
            end
            
            
        end
    end
end
